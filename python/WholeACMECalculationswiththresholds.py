#-------------------------------------------------------------------------------
# Name:        WholeACMEcalculationsWithThresholds
# Purpose:     Derive all cirque metrics based on Cirque outlines, a DEM, and user-provided thresholds
#
# Author: Yingkui Li
# This program derive cirque related metrics based on cirque outlines and a DEM
# The first step is to determine the cirque threshold points
# The second step is to derive length and width info, as well as the area, and parameters
# The third step is to detive the 3D statistics and hypsometric parameters
# Some of the codes are revised based on the ACME codes by Ramon Pellitero and Matteo Spagnolo 2016
# 
# Created:     05/26/2023
# Copyright:   (c) Yingkui Li 2023
#-------------------------------------------------------------------------------

from __future__ import division
import arcpy
from arcpy import env
from arcpy.sa import *
import math
import time
import numpy as np
from scipy.optimize import curve_fit
from scipy import optimize

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"

arcpy.Delete_management("in_memory") ### Empty the in_memory
ArcGISPro = 0
arcpy.AddMessage("The current python version is: " + str(sys.version_info[0]))
if sys.version_info[0] == 2:  ##For ArcGIS 10, need to check the 3D and Spatial Extensions
    try:
        if arcpy.CheckExtension("Spatial")=="Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"

    try:
        if arcpy.CheckExtension("3D")=="Available":
            arcpy.CheckOutExtension("3D")
        else:
            raise Exception ("not extension available")
            #print "not extension available"
    except:
        raise Exception ("unable to check out extension")
        #print "unable to check out extension"
elif sys.version_info[0] == 3:  ##For ArcGIS Pro
    ArcGISPro = 1
    #pass ##No need to Check
else:
    raise Exception("Must be using Python 2.x or 3.x")
    exit()   


# Polynomial Regression
def polyfit(x, y, degree):
    results = {}

    coeffs = np.polyfit(x, y, degree)

     # Polynomial Coefficients
    results['polynomial'] = coeffs.tolist()

    # r-squared
    p = np.poly1d(coeffs)
    # fit values, and mean
    yhat = p(x)                         # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    results['determination'] = ssreg / sstot

    return results

def k_curve (x, c):
    return (1-x) * np.exp(c * x)

# normalized K_curve Regression
def k_curve_fit(x, y):
    popt, pcov = curve_fit(k_curve, x, y)
    c = popt[0]
    # r-squared
    yhat = k_curve(x, c)             # or [p(z) for z in x]
    ybar = np.sum(y)/len(y)          # or sum(y)/len(y)
    ssreg = np.sum((yhat-ybar)**2)   # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar)**2)    # or sum([ (yi - ybar)**2 for yi in y])
    R2 = ssreg / sstot
    if R2 > 1:
        R2 = 1/R2
    return (c, R2)

#---------------------------------------------------------------------------------------------------------------
# This function calculates the distance between two points
#--------------------------------------------------------------------------------------------------------------- 
def Dist(x1,y1,x2,y2):
    return math.sqrt(math.pow(math.fabs(x1-x2),2)+math.pow(math.fabs(y1-y2),2))

#------------------------------------------------------------------------------------------------------------
# This function check each line in the line feature and make sure the line is from low elevation to high
# elevation (Glacier flowline needs from low to high elevation in order to reconstruct the paleo ice thickness).
# It is revised from the codes by Pellitero et al.(2016) in GlaRe.
#------------------------------------------------------------------------------------------------------------
def Check_If_Flip_Line_Direction(line, dem):
    arcpy.AddField_management(line, "Flip", "Long", "", "", "", "", "", "", "")
    updCursor = arcpy.da.UpdateCursor(line, ["Shape@", "Flip"])
    nFlip = 0
    i = 0
    for curLine in updCursor:
        Startpoint2 = curLine[0].firstPoint
        Startpoint = arcpy.Geometry('Point',Startpoint2)
        coord= str(Startpoint2.X)+" "+str(Startpoint2.Y)
        Cellvalue=arcpy.GetCellValue_management(dem, coord)
        Startpoint.Z=Cellvalue.getOutput(0)
        Lastpoint2 = curLine[0].lastPoint
        Lastpoint = arcpy.Geometry('Point',Lastpoint2)
        coord= str(Lastpoint2.X)+" "+str(Lastpoint2.Y)
        Cellvalue=arcpy.GetCellValue_management(dem, coord)
        Lastpoint.Z=Cellvalue.getOutput(0)
        if Startpoint.Z >= Lastpoint.Z:  ##Flip = True use equal in case the start and end point are the same
            nFlip = nFlip + 1
            curLine[1] = 1
        else:  ##Flip = False
            curLine[1] = 0
        updCursor.updateRow(curLine)
        i += 1 
        
    if nFlip > 0:
        arcpy.MakeFeatureLayer_management(line, "lyrLines")
        arcpy.SelectLayerByAttribute_management("lyrLines", "NEW_SELECTION", '"Flip" > 0')
        #arcpy.AddMessage("Count of selected features is " + str(nFlip))
        arcpy.FlipLine_edit("lyrLines")  ##Need to change to lyrLines
        arcpy.SelectLayerByAttribute_management("lyrLines", "CLEAR_SELECTION")

    arcpy.DeleteField_management (line, "Flip")
    
    del updCursor
    if i>0:
        del curLine

 
def getAngle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang + 360 if ang < 0 else ang

def orientation(p1, p2, p3):
     
    # to find the orientation of
    # an ordered triplet (p1,p2,p3)
    # function returns the following values:
    # 0 : Collinear points
    # 1 : Clockwise points
    # 2 : Counterclockwise
    val = (float(p2[1] - p1[1]) * (p3[0] - p2[0])) - float((p2[0] - p1[0]) * (p3[1] - p2[1]))
    if (val > 0):
        # Clockwise orientation
        return 1
    elif (val < 0):
        # Counterclockwise orientation
        return 2
    else:
        # Collinear orientation
        return 0

def line_direction(startpnt, endpnt): ##Calculate the line direction clockwise from north (0-360)
    dx  = endpnt[0] - startpnt[0]
    dy  = endpnt[1] - startpnt[1]

    angle = 180.0/math.pi * math.atan2(dy, dx)
    angle = 90 - angle
    if angle < 0:
        angle += 360
    return angle

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]
    try:
        div = det(xdiff, ydiff)
    except:
        div = 0
        
    if div == 0:
        #raise Exception('lines do not intersect')
        #print ('lines do not intersect')
        return 0,(0,0)
    

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return 1, (x, y)

def connect_line_sections(line):
    lineArray = arcpy.da.FeatureClassToNumPyArray(line,"SHAPE@LENGTH")
    if len(lineArray) == 1:
        return line
    startpnts = "in_memory\\startpnts"
    endpnts = "in_memory\\endpnts"
    arcpy.FeatureVerticesToPoints_management(line, startpnts, "START")
    pntArray = arcpy.da.FeatureClassToNumPyArray(startpnts,["SHAPE@X", "SHAPE@Y"])
    startpntX = np.array([item[0] for item in pntArray])
    startpntY = np.array([item[1] for item in pntArray])
    
    arcpy.FeatureVerticesToPoints_management(line, endpnts, "End")
    pntArray = arcpy.da.FeatureClassToNumPyArray(endpnts,["SHAPE@X", "SHAPE@Y"])
    endpntX = np.array([item[0] for item in pntArray])
    endpntY = np.array([item[1] for item in pntArray])

    distlist = []
    pntCount = len(lineArray)
    for i in range(pntCount):
        for j in range(pntCount):
            if i != j:
                distance = Dist (startpntX[i], startpntY[i], endpntX[j], endpntY[j])
            else:
                distance = 99999

            distlist.append(distance)

    #arcpy.AddMessage(distlist)
    ##Find the shortest distance inorder until pntCount - 1
    ##sort the distlist from small to large
    distArr = np.array(distlist)
    indexArr = np.argsort(distArr)
    #arcpy.AddMessage(indexArr)

    connectionline = arcpy.CreateFeatureclass_management("in_memory", "connectionline","POLYLINE","","","",line)
    new_line_cursor = arcpy.da.InsertCursor(connectionline, ('SHAPE@'))
    
    for i in range(pntCount-1):
        index = indexArr[i]
        ##Figure out the start and end location of the points
        start = int(index / pntCount)
        end = index % pntCount
        array = arcpy.Array([arcpy.Point(startpntX[start],startpntY[start]),arcpy.Point(endpntX[end], endpntY[end])])
        polyline = arcpy.Polyline(array)
        new_line_cursor.insertRow([polyline])

    del new_line_cursor

    arcpy.Append_management(connectionline, line, "NO_TEST")
    arcpy.Dissolve_management(line, "in_memory\\line_dissolve")
    arcpy.CopyFeatures_management("in_memory\\line_dissolve", line)

    arcpy.Delete_management ("in_memory\\line_dissolve")
    arcpy.Delete_management (connectionline)
    
    return line

def define_circle_center(p1, p2, p3):
    """
    Returns the center and radius of the circle passing the given 3 points.
    In case the 3 points form a line, returns the middle point.
    """
    temp = p2[0] * p2[0] + p2[1] * p2[1]
    bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
    cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
    det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])

    if abs(det) < 1.0e-6:
        #arcpy.AddMessage("The three points are along the same line")
        cx = (p1[0] + p2[0]) / 2
        cy = (p1[1] + p2[1]) / 2
        return (cx, cy)
        #return None
    
    # Center of circle
    cx = (bc*(p2[1] - p3[1]) - cd*(p1[1] - p2[1])) / det
    cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det
    
    #radius = np.sqrt((cx - p1[0])**2 + (cy - p1[1])**2)

    return (cx, cy)

#------------------------------------------------------------------------------------------------------------
# This function calculates the plan closure for a cirque. The codes are revised basedd on ACME codes
# This the revised new codes for the plan closure of the original ACME
#------------------------------------------------------------------------------------------------------------
def plan_closSPA(cirqueDEM):
    #tool to get the plan_closure
    mid_height = (cirqueDEM.maximum + cirqueDEM.minimum) / 2
    #arcpy.AddMessage("mid_height is: " + str(mid_height))
    
    #Step 2: Get the contourline of the mid-height
    midHcontour = Contour(cirqueDEM, "in_memory/cont", 10000, mid_height)
    lineArray = arcpy.da.FeatureClassToNumPyArray("in_memory/cont","SHAPE@LENGTH")
    if len(lineArray) > 1: ##if more than one contour lines
        arcpy.AddMessage("multiple contour sections, connect into one!")
        connect_line_sections("in_memory/cont")
    
    lineArray = arcpy.da.FeatureClassToNumPyArray("in_memory/cont","SHAPE@LENGTH")
    lengthArr = np.array([item[0] for item in lineArray])
    max_length = np.max(lengthArr)
    with arcpy.da.UpdateCursor("in_memory/cont", "SHAPE@LENGTH") as cursor:
        for row in cursor:
            if row[0] < (max_length-0.1): ##This step reomve the small contour lines
                arcpy.AddMessage("Delete small contour lines")
                cursor.deleteRow()
    del cursor, row

    #Step 3: find the coordinates of first, mid and last points along the mid_height contour
    ##Get the start and end points of the contour line
    with arcpy.da.SearchCursor("in_memory/cont", "SHAPE@") as cursor:
        for row in cursor:
            #Get start point
            startpt = row[0].firstPoint
            x_start = startpt.X
            y_start = startpt.Y

            #Get end point
            endpt = row[0].lastPoint
            x_end = endpt.X
            y_end = endpt.Y    

    #Get the center point of contourline
    arcpy.FeatureVerticesToPoints_management("in_memory/cont", "in_memory/mid_point", "MID")
    pntArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\mid_point",["SHAPE@X", "SHAPE@Y"])
    pntXArr = np.array([item[0] for item in pntArray])
    pntYArr = np.array([item[1] for item in pntArray])
    x_mid = np.mean(pntXArr)
    y_mid = np.mean(pntYArr)
    
    #Step 4: find center point of the three points by circle
    centre = define_circle_center((x_start, y_start), (x_end, y_end), (x_mid, y_mid))
    
    #define end_start and mid_centre segments as polylines
    array_centre_mid=arcpy.Array([arcpy.Point(x_mid, y_mid),arcpy.Point(centre[0], centre[1])])
    segment_centre_mid=arcpy.Polyline(array_centre_mid)
    array_start_end=arcpy.Array([arcpy.Point(x_start, y_start),arcpy.Point(x_end, y_end)])
    segment_start_end=arcpy.Polyline(array_start_end)

    #Step 5: Find the angle of the start point, centre point and the end point
    Angle = getAngle((x_start, y_start), centre , (x_end, y_end))
    #arcpy.AddMessage(str(Angle))
    #verify whether the mid_centre segment intersect end_start segment
    if (segment_centre_mid.crosses(segment_start_end)==False):
        #calculate 360 degrees - the angle between centre, start and end points
        #Angle = ((2*math.pi - (math.acos(((math.pow(dist_centre_end,2)+math.pow(dist_centre_start,2)-math.pow(dist_start_end,2))/(2*dist_centre_end*dist_centre_start)))))*180/math.pi)
        if Angle < 180:
            Angle = 360 - Angle
    else:
        #calculate the angle between centre, start and end points
        #Angle = (((math.acos(((math.pow(dist_centre_end,2)+math.pow(dist_centre_start,2)-math.pow(dist_start_end,2))/(2*dist_centre_end*dist_centre_start)))))*180/math.pi)
        if Angle > 180:
            Angle = 360 - Angle
                
    ##delete the temp dataset
    arcpy.Delete_management ("in_memory/cont")
    arcpy.Delete_management ("in_memory/mid_point")

    return Angle   

'''
##old codes from ACME
def plan_closSPA(cirqueDEM):
    #tool to get the plan_closure
    meanH = int(cirqueDEM.mean)
    ##print meanH
    mid_height=(meanH)
    midHcontour = Contour(cirqueDEM, "in_memory/cont", 10000, mid_height)
    geometry = arcpy.CopyFeatures_management(midHcontour, arcpy.Geometry())
    #get total length of geometry to select only the longest contour if at that elevation there are more than one line
    lenghtlist=[]
    for x in geometry:
        #arcpy.AddMessage("length is:" + str(x.length))
        if x.length > 100: ##contour needs to > 100 m
            lenghtlist.append(x.length)
    if len(lenghtlist) > 0:
        maxlength=max(lenghtlist)
        index=lenghtlist.index(maxlength)
        goodgeometry=geometry[index]
        maximum=int(geometry[index].length)
        #find the coordinates of first, mid and last points along the mid_height contour
        point_coord_list=[]
        for m in range(0, maximum+1, int(maximum/2)):
            points = goodgeometry.positionAlongLine (m)
            centroid = points.centroid
            point_coord_list.append(centroid.X)
            point_coord_list.append(centroid.Y)
        #define the coordinates
        x_start,y_start=point_coord_list[0],point_coord_list[1]
        x_end,y_end=point_coord_list[4],point_coord_list[5]
        x_mid,y_mid=point_coord_list[2],point_coord_list[3]

        #to avoid dividing by 0 in the next section, i.e. when getting s1 and s2, x_mid and x_end and x1 cannot have the same value
        if x_end==x_start:
            x_end+=0.00001
        elif x_mid==x_end:
            x_end+=0.00001
        else:
            pass

    
        #end_start and mid_end lines midpoints coordinates

        m1x,m1y= (x_end+x_start)/2, (y_end+y_start)/2
        m2x,m2y= (x_mid+x_end)/2, (y_mid+y_end)/2

        # slope of the end_start and mid_end lines
        s1=(y_end-y_start)/(x_end-x_start)
        if s1 == 0:
            s1 = 0.00001
        ##print s1
        s2=(y_mid-y_end)/(x_mid-x_end)
        ##print s2
        if s2 == 0:
            s2 = 0.00001

        #inverse slope
        is1=-1*(1/s1)
        is2=-1*(1/s2)

        #equations that enable to define the point of intersection between the two lines passing by the midpoints and perpendicular to the end_start and mid_end segments
        a=np.array([[-is1,1],[-is2,1]])
        b=np.array([(-m1x*is1+m1y),(-m2x*is2+m2y)])
        try:
            centre=np.linalg.solve(a,b)
        except:
            arcpy.AddMessage("three points are probably colinear")
            #print "three points are probably colinear"
            return 0

        #measure distances between key points
        dist_centre_start = math.sqrt((math.pow((x_start-centre[0]),2))+(math.pow((y_start-centre[1]),2)))
        dist_centre_end = math.sqrt((math.pow((x_end-centre[0]),2))+(math.pow((y_end-centre[1]),2)))
        dist_centre_mid = math.sqrt((math.pow((x_mid-centre[0]),2))+(math.pow((y_mid-centre[1]),2)))
        dist_start_end = math.sqrt((math.pow((x_start-x_end),2))+(math.pow((y_start-y_end),2)))
        #define end_start and mid_centre segments as polylines
        array_centre_mid=arcpy.Array([arcpy.Point(x_mid, y_mid),arcpy.Point(centre[0], centre[1])])
        segment_centre_mid=arcpy.Polyline(array_centre_mid)
        array_start_end=arcpy.Array([arcpy.Point(x_start, y_start),arcpy.Point(x_end, y_end)])
        segment_start_end=arcpy.Polyline(array_start_end)
        #verify whether the mid_centre segment intersect end_start segment
        if segment_centre_mid.crosses(segment_start_end)==False:
            #calculate 360 degrees - the angle between centre, start and end points
            Angle = ((2*math.pi - (math.acos(((math.pow(dist_centre_end,2)+math.pow(dist_centre_start,2)-math.pow(dist_start_end,2))/(2*dist_centre_end*dist_centre_start)))))*180/math.pi)
        else:
            #calculate the angle between centre, start and end points
            Angle = (((math.acos(((math.pow(dist_centre_end,2)+math.pow(dist_centre_start,2)-math.pow(dist_start_end,2))/(2*dist_centre_end*dist_centre_start)))))*180/math.pi)
    else:
        Angle = 0

    ##delete the temp dataset
    arcpy.Delete_management ("in_memory/cont")

    return Angle
'''

#------------------------------------------------------------------------------------------------------------
# This function calculates the plan closure for a cirque based on the midAlt contours. The plan closure is determined as the
# difference between the line directions of the start and end 100 m sections of the contour line
# This method is based on IS Evans' method
#------------------------------------------------------------------------------------------------------------
def plan_closISE(cirqueDEM, midaltcontur): ##, startline, endline):
    #tool to get the plan_closure
    ##Step 1: get the mid height by mean elevation or by mid of maximum or minimum??
    Z_min = int(cirqueDTM.minimum)
    Z_max = int(cirqueDTM.maximum)
    mid_height = (Z_min + Z_max) / 2

    #arcpy.AddMessage("mid_alt is: " + str(mid_height))

    #Step 2: Get the contourline of the mid-height
    cont = "in_memory\\cont"
    midHcontour = Contour(cirqueDEM, cont, 10000, float(mid_height))

    lineArray = arcpy.da.FeatureClassToNumPyArray(cont,"SHAPE@LENGTH")
    if len(lineArray) > 1: ##if more than one contour lines
        arcpy.AddMessage("multiple contour sections, connect into one!")
        connect_line_sections(cont)
    
    ##Copy the contour out
    arcpy.CopyFeatures_management(cont, midaltcontur)

    #Step 3: define the perpendicular line for the start point of the mid_height contour (100 m buffer intersection line with the midalt contour)
    startpnt = "in_memory\\startpnt"
    bothpnts = "in_memory\\bothpnts"
    arcpy.FeatureVerticesToPoints_management(cont, startpnt, "START")
    
    tmpbuf = "in_memory\\tmpbuf"
    arcpy.Buffer_analysis(startpnt, tmpbuf, "100 Meters")
    cont_clip = "in_memory\\cont_clip"
    arcpy.Clip_analysis(cont, tmpbuf, cont_clip)
    clip_cont_single = "in_memory\\clip_cont_single"
    arcpy.MultipartToSinglepart_management(cont_clip, clip_cont_single)
    lineArray = arcpy.da.FeatureClassToNumPyArray(clip_cont_single,"SHAPE@LENGTH")

    if len(lineArray) > 1:
        lengthArr = np.array([item[0] for item in lineArray])
        #arcpy.AddMessage(lengthArr)
        
        max_length = np.max(lengthArr)
        with arcpy.da.UpdateCursor(clip_cont_single, "SHAPE@LENGTH") as cursor:
            for row in cursor:
                if row[0] < (max_length-0.1): ##This step reomve the small contour lines
                    arcpy.AddMessage("Delete small contour lines")
                    cursor.deleteRow()
        del cursor, row

    arcpy.FeatureVerticesToPoints_management(clip_cont_single, bothpnts, "BOTH_ENDS")
    pntArray = arcpy.da.FeatureClassToNumPyArray(bothpnts,["SHAPE@X", "SHAPE@Y"])
    pntXArr = np.array([item[0] for item in pntArray])
    pntYArr = np.array([item[1] for item in pntArray])

    startx = pntXArr[0]
    starty = pntYArr[0]
    startpnt2x = pntXArr[1]
    srartpnt2y = pntYArr[1]

    start_angle = line_direction((startx,starty), (startpnt2x, srartpnt2y))
    #arcpy.AddMessage("Start_angle is: " + str(start_angle))

    #startperpline = "in_memory\\startperpline"
    #extend_line((startx,starty), (startpnt2x, srartpnt2y), 20000, cont, startperpline)
    ##Copy startline
    #arcpy.CopyFeatures_management(startperpline, startline)

    #Step 4: define the perpendicular line for the end point of the mid_height contour (100 m buffer intersection line with the midalt contour)
    endpnt = "in_memory\\endpnt"
    arcpy.FeatureVerticesToPoints_management(cont, endpnt, "END")

    arcpy.Buffer_analysis(endpnt, tmpbuf, "100 Meters")
    arcpy.Clip_analysis(cont, tmpbuf, cont_clip)
    arcpy.MultipartToSinglepart_management(cont_clip, clip_cont_single)
    lineArray = arcpy.da.FeatureClassToNumPyArray(clip_cont_single,"SHAPE@LENGTH")

    if len(lineArray) > 1:
        lengthArr = np.array([item[0] for item in lineArray])
        #arcpy.AddMessage(lengthArr)
        max_length = np.max(lengthArr)
        with arcpy.da.UpdateCursor(clip_cont_single, "SHAPE@LENGTH") as cursor:
            for row in cursor:
                if row[0] < (max_length-0.1): ##This step reomve the small contour lines
                    arcpy.AddMessage("Delete small contour lines")
                    cursor.deleteRow()
        del cursor, row

    arcpy.FeatureVerticesToPoints_management(clip_cont_single, bothpnts, "BOTH_ENDS")
    pntArray = arcpy.da.FeatureClassToNumPyArray(bothpnts,["SHAPE@X", "SHAPE@Y"])
    pntXArr = np.array([item[0] for item in pntArray])
    pntYArr = np.array([item[1] for item in pntArray])

    endpnt2x = pntXArr[0]
    endpnt2y = pntYArr[0]
    endx = pntXArr[1]
    endy = pntYArr[1]

    end_angle = line_direction((endpnt2x, endpnt2y), (endx, endy))
    #arcpy.AddMessage("end_angle is: " + str(end_angle))


    #endperpline = "in_memory\\endperpline"
    #extend_line((endpnt2x, endpnt2y), (endx, endy), 20000, cont, endperpline)
    ##Copy endline
    #arcpy.CopyFeatures_management(endperpline, endline)

    #Step 5: Get the intersection point between the two perplines
    intersect_result = line_intersection(((startx,starty), (startpnt2x, srartpnt2y)), ((endpnt2x, endpnt2y), (endx, endy)))
    
    angleType = 0
    bFlag = 0
    Angle = 0
    #if len(pntXArr) > 0:
    if intersect_result[0] > 0:
        intersectx = intersect_result[1][0]
        intersecty = intersect_result[1][1]

        ##determine the distance between start line section to the intersect 
        start_dist_diff = Dist(startx,starty,intersectx,intersecty) - Dist(startpnt2x,srartpnt2y,intersectx,intersecty)
        ##determine the distance between end line section to the intersect 
        end_dist_diff = Dist(endx,endy,intersectx,intersecty) - Dist(endpnt2x,endpnt2y,intersectx,intersecty)

        flag = start_dist_diff * end_dist_diff
        if flag > 0:
            bFlag = 0
        else:
            bFlag = 1
        
        #Step 5: Find the angle of the start point, centre point and the end point
        angleType = orientation((startx, starty), (intersectx, intersecty) , (endx, endy))
        Angle = end_angle - start_angle
        if Angle < 0:
            Angle += 360
        #arcpy.AddMessage("Angle type is: " + str(angleType))
        #arcpy.AddMessage("Flag is: " + str(bFlag))

        if (angleType == 1):
            if (flag > 0):
                if start_dist_diff < 0:
                    Angle = 360 - Angle 
            else:
                if Angle > 180:
                    Angle = 360 - Angle 
        
        if (angleType == 2):
            if (flag > 0):
                if start_dist_diff > 0:
                    Angle = 360 - Angle 
            else:
                if Angle > 180:
                    Angle = 360 - Angle 
                
        #arcpy.AddMessage("Plan closure is: " + str(Angle))
    else:
        Angle = 180

    ##delete the temp dataset
    arcpy.Delete_management (cont)
    arcpy.Delete_management (startpnt)
    arcpy.Delete_management (endpnt)
    arcpy.Delete_management (bothpnts)
    arcpy.Delete_management (tmpbuf)
    arcpy.Delete_management (cont_clip)
    arcpy.Delete_management (clip_cont_single)
    #arcpy.Delete_management (startperpline)
    #arcpy.Delete_management (endperpline)

    return (Angle, start_angle, end_angle, angleType, bFlag)

##Main program
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)
InputThresholds = arcpy.GetParameterAsText(2)
OutputLength = arcpy.GetParameterAsText(3)
OutputWidth = arcpy.GetParameterAsText(4)
OutputMidAltContour = arcpy.GetParameterAsText(5)
outputCirques = arcpy.GetParameterAsText(6)

#environments

spatialref=arcpy.Describe(InputCirques).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

arcpy.Delete_management("in_memory") ### Empty the in_memory

##Obtain the cellsize of the DEM
dem = arcpy.Raster(InputDEM)
cellsize = dem.meanCellWidth
cellsize_int = int(cellsize)

cirques_copy = arcpy.env.scratchGDB + "\\cirques_copy"
arcpy.CopyFeatures_management(InputCirques, cirques_copy)

#a "list" where the name of fields from the attributed table are copied in
Fieldlist=[]
ListFields=arcpy.ListFields(cirques_copy)
for x in ListFields:
    Fieldlist.append(x.baseName)

##Add new fields
if "ID_cirque" in Fieldlist:
    pass
else:
    arcpy.AddField_management(cirques_copy, "ID_cirque", 'Long')

arcpy.CalculateField_management(cirques_copy,"ID_cirque",str("!"+str(arcpy.Describe(cirques_copy).OIDFieldName)+"!"),"PYTHON_9.3")

##Map projection and DEM resolution info
if "Projection" in Fieldlist:
    pass
else:
    arcpy.AddField_management(cirques_copy, "Projection", 'Text')

if "DEMresolution" in Fieldlist:
    pass
else:
    arcpy.AddField_management(cirques_copy, "DEMresolution", "DOUBLE")

if "FocusMethod" in Fieldlist:
    pass
else:
    arcpy.AddField_management(cirques_copy, "FocusMethod", 'Text')

##location variables
new_fields = ("Lon","Lat","Easting","Northing") ##All double variables count = 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 4)


##Size variables
new_fields = ("L","W","H","CS","Perimeter", "A2D", "A3D") ##All Integer variables count = 7
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "LONG",10)

##Shape variables
new_fields = ("L_W","L_H","W_H","A3D_A2D","Circular") ##All float variables count = 5
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 3)

new_fields = ("Slope_mean", "Slope_max", "Slope_min", "Slpgt33","Slplt20","Slp20to33", "Aspectmean") ##All float variables count = 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 1)

new_fields = ("Asp_east","Asp_north")##All float variables count = 2
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 3)

new_fields = ("Plan_closSPA", "Plan_closISE", "Prof_clos") ##All float variables count = 3
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 1)

##Altitude variables
new_fields = ("Z_min","Z_max", "Z_mean","Z_median","Z_mid", "Hypsomax") ##All Integer variables count = 6
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "LONG",10)

if "HI" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirques_copy, "HI", "DOUBLE", 10,3)

##Axis variables
new_fields = ("Axprofclos", "Axhli", "Axasp","Axgrad") ## All float variables 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 2)

if "Axamp" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirques_copy, "Axamp", "LONG", 10)


##axis curve-fit variables 
new_fields = ("L_Exp_A","L_Exp_B","L_Exp_R2","L_Kcurv_C","L_Kcurv_R2","W_Quad_C", "W_Quad_R2") ##float variables with high digits
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "FLOAT",10, 4)

##Catchment variables
if "Maxabalt" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirques_copy, "Maxabalt", "LONG", 10)

if "Pctabarea" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirques_copy, "Pctabarea", "DOUBLE",10, 2)

##Step 1: derive Location variables
arcpy.AddMessage("Step 1: Deriving Location variables: Lat, Long, Northing (km) and Easting (km)")
'''
new_fields = ("M_East","M_North")
##Get the easting and northing
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 4)

spatial_ref = arcpy.Describe(cirques_copy).spatialReference
#arcpy.AddMessage(spatial_ref.name)
#arcpy.AddMessage(spatial_ref.XYResolution)
if "GCS" in spatial_ref.name:
    ##Need to convert the map projection to a UTM or similar
    arcpy.AddMessage("The projection is GCS. Need to reproject the file to a UTM or other projected coordinate system")
    exit()
arcpy.management.CalculateGeometryAttributes(cirques_copy, [["M_East", "INSIDE_X"],["M_North", "INSIDE_Y"]])

fields = ("M_East", "Easting", "M_North", "Northing", "Projection", "DEMresolution", "FocusMethod")
method = "User-specified"
with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
    for row in cursor:
       row[1] = row[0]/1000 ##Convert to km
       row[3] = row[2]/1000 ##Convert to km
       row[4] = spatial_ref.name
       row[5] = cellsize
       row[6] = method
       cursor.updateRow(row)
del row, cursor
arcpy.DeleteField_management(cirques_copy,["M_East", "M_North"])

##get the lat and long
#wkt = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],\
#              PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];\
#              -400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;\
#              0.001;0.001;IsHighPrecision"
#sr = arcpy.SpatialReference(text=wkt)

sr = arcpy.SpatialReference(4326) ##WGS 84
arcpy.management.CalculateGeometryAttributes(cirques_copy, [["Lon", "INSIDE_X"],["Lat", "INSIDE_Y"]], "", "", sr, "DD")
'''
spatial_ref = arcpy.Describe(cirques_copy).spatialReference
#arcpy.AddMessage(spatial_ref.name)
#arcpy.AddMessage(spatial_ref.XYResolution)
if "GCS" in spatial_ref.name:
    ##Need to convert the map projection to a UTM or similar
    arcpy.AddMessage("The projection is GCS. Need to reproject the file to a UTM or other projected coordinate system")
    exit()

##Obtain the cellsize of the DEM
dem = arcpy.Raster(InputDEM)
cellsize = dem.meanCellWidth
#arcpy.AddMessage("Cell size: " + str(cellsize))

if ArcGISPro:
    new_fields = ("M_East","M_North")
    ##Get the easting and northing
    for field in new_fields:
        if field in Fieldlist:
            pass
        else:
            arcpy.AddField_management(cirques_copy, field, "DOUBLE",10, 4)

    arcpy.management.CalculateGeometryAttributes(cirques_copy, [["M_East", "INSIDE_X"],["M_North", "INSIDE_Y"]])

    fields = ("M_East", "Easting", "M_North", "Northing", "Projection", "DEMresolution", "FocusMethod")
    method = "User-specified"
    with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
        for row in cursor:
           row[1] = row[0]/1000 ##Convert to km
           row[3] = row[2]/1000 ##Convert to km
           row[4] = spatial_ref.name
           row[5] = cellsize
           row[6] = method
           cursor.updateRow(row)
    del row, cursor
    arcpy.DeleteField_management(cirques_copy,["M_East", "M_North"])

else: #For ArcGIS 10
    arcpy.AddGeometryAttributes_management(cirques_copy, "CENTROID_INSIDE")
    fields = ("INSIDE_X", "Easting", "INSIDE_Y", "Northing", "Projection", "DEMresolution", "FocusMethod")
    method = "User-specified"
    with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
        for row in cursor:
           row[1] = row[0]/1000 ##Convert to km
           row[3] = row[2]/1000 ##Convert to km
           row[4] = spatial_ref.name
           row[5] = cellsize
           row[6] = method
           cursor.updateRow(row)
    del row, cursor
    arcpy.DeleteField_management(cirques_copy,["INSIDE_X", "INSIDE_Y"])

##get the lat and long
#wkt = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],\
#              PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];\
#              -400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119522E-09;\
#              0.001;0.001;IsHighPrecision"
#sr = arcpy.SpatialReference(text=wkt)

sr = arcpy.SpatialReference(4326) ##WGS 84
if ArcGISPro:
    arcpy.management.CalculateGeometryAttributes(cirques_copy, [["Lon", "INSIDE_X"],["Lat", "INSIDE_Y"]], "", "", sr, "DD")
else: ##For ArcGIS 10
    arcpy.AddGeometryAttributes_management(cirques_copy, "CENTROID_INSIDE", "", "", sr)
    fields = ("INSIDE_X", "Lon", "INSIDE_Y", "Lat")
    with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
        for row in cursor:
           row[1] = row[0]
           row[3] = row[2]
           cursor.updateRow(row)
    del row, cursor
    arcpy.DeleteField_management(cirques_copy,["INSIDE_X", "INSIDE_Y"])


    Field_list = []
    List_Fields = arcpy.ListFields(cirques_copy)
    for x in List_Fields:
        Field_list.append(x.baseName)

    ##Delete other inside variables created
    if "INSIDE_Z" in Field_list:
        arcpy.DeleteField_management(cirques_copy,["INSIDE_Z"])
    if "INSIDE_M" in Field_list:
        arcpy.DeleteField_management(cirques_copy,["INSIDE_M"])

##Step 2: derive Perimeter, Area_2D, and Circularity
arcpy.AddMessage("Step 2: Deriving Perimeter, Area_2D, and Circularity")

fields = ("Perimeter", "A2D", "Circular", "SHAPE@LENGTH", "SHAPE@AREA")
with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
    for row in cursor:
       row[0]=row[3]
       row[1]=row[4]
       row[2]=row[3]/(2*math.pi*((math.sqrt(row[4]/math.pi))))
       cursor.updateRow(row)
del row, cursor

##Step 3: derive Length and width-related metrics 
arcpy.AddMessage("Step 3: Deriving length and width-related metrics")

cirque_length = arcpy.CreateFeatureclass_management("in_memory", "cirque_length", "POLYLINE", "","","",spatialref)
#add ID field to attribute table of new shapefile to be matched with the FID of the original cirque polygon shapefile
arcpy.AddField_management(cirque_length, "ID_cirque", "INTEGER", 6)
cirque_width = arcpy.CreateFeatureclass_management("in_memory", "cirque_width", "POLYLINE", "","","",spatialref)
#add ID field to attribute table of new shapefile to be matched with the FID of the original cirque polygon shapefile
arcpy.AddField_management(cirque_width, "ID_cirque", "INTEGER", 6)

#routine for all cirques and threshold thresh_points
#CirqueID = arcpy.Describe(cirques_copy).OIDFieldName
CirqueID = "ID_cirque"
with arcpy.da.UpdateCursor(cirques_copy, ["L", "W", "L_W", "SHAPE@", CirqueID]) as cursor:
    for row in cursor:

        arcpy.CopyFeatures_management(row[3], "in_memory/feature")
        ##Need to use the outline to determine the cirque threshold points
        arcpy.PolygonToLine_management("in_memory/feature", "in_memory\\feature_lines")
        thresh_p = arcpy.SpatialJoin_analysis(InputThresholds, "in_memory\\feature_lines", "in_memory\\thres", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "30 Meters", "#")
        n_thresholds = arcpy.GetCount_management(thresh_p)#checking that there is only one intersection
        count=int(n_thresholds.getOutput(0)) #checking that there is only one intersection
        i = 2
        while count < 1:
            ##increase the search distance
            search_dis = str(i*30) + " Meters"
            #arcpy.AddMessage(search_dis)
            thresh_p = arcpy.SpatialJoin_analysis(InputThresholds, "in_memory\\feature_lines", "in_memory\\thres", "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", search_dis, "#")
            n_thresholds = arcpy.GetCount_management(thresh_p)#checking that there is only one intersection
            count=int(n_thresholds.getOutput(0)) #checking that there is only one intersection
            i += 1

        cirque_id = row[4]
        arcpy.AddMessage("Processing cirque #" + str(cirque_id))


        if count > 1: ##if there are multiple thresholds, take the lowest one
            #arcpy.AddMessage("The number of threshold points: " + str(count))
            ExtractValuesToPoints(thresh_p, InputDEM, "in_memory\\thres_with_Ele", "INTERPOLATE", "VALUE_ONLY") 
            pntArr = arcpy.da.FeatureClassToNumPyArray("in_memory\\thres_with_Ele", "RASTERVALU")
            elevs = np.array([item[0] for item in pntArr])
            min_elev = min(elevs)

            with arcpy.da.UpdateCursor("in_memory\\thres_with_Ele", "RASTERVALU") as cursor2:
                for row2 in cursor2:
                    if row2[0] > min_elev: ##This step reomve the small contour lines
                        #arcpy.AddMessage("Delete one threshold lines")
                        cursor2.deleteRow()
            del cursor2, row2
            arcpy.CopyFeatures_management("in_memory\\thres_with_Ele", thresh_p)
            
        
        #if count==1:
        cirgeom=arcpy.CopyFeatures_management(row[3], arcpy.Geometry())#geometry for the cirque
        centroid=row[3].centroid#cirque centroid (needs ArcGIS 10.1)
        thresh_p_geometry = arcpy.CopyFeatures_management(thresh_p, arcpy.Geometry())
        firstpoint=thresh_p_geometry[0].firstPoint
        V1X=(firstpoint.X-centroid.X)+0.0001 #first coordinate of vector center-first point
        V1Y=(firstpoint.Y-centroid.Y)+0.0001 #second coordinate of vector center-first point
        slope=V1Y/V1X #line slope
        distance = math.sqrt(math.pow(V1X,2)+math.pow(V1Y,2))
        pointxlength=(distance) * 10
        pointylength=(slope)*(distance) * 10
        pointXpos=centroid.X+pointxlength
        pointYpos=centroid.Y+pointylength
        pointA=arcpy.Point(pointXpos, pointYpos)#first corner of cross line
        pointXneg=centroid.X-pointxlength
        pointYneg=centroid.Y-pointylength
        pointB=arcpy.Point(pointXneg, pointYneg) #second corner of cross line
        lengtharray1=arcpy.Array([pointA,centroid]) #first cross line
        LENGTH1=arcpy.Geometry("Polyline",lengtharray1, spatialref) #first cross line geometry
        lengtharray2=arcpy.Array([pointB,centroid]) #second cross line
        LENGTH2=arcpy.Geometry("Polyline",lengtharray2, spatialref) #second cross line geometry

        #turn boundary into geometry
        bound=cirgeom[0].boundary()#cirque boundary (needs ArcGIS 10.1)
        array=bound.getPart()
        boundary=arcpy.Geometry("Polyline",array,spatialref)

        point1=LENGTH1.intersect(boundary,1) # first intersection of boundary and cross line
        point2=LENGTH2.intersect(boundary,1) # second intersection of boundary and cross line

        #length line
        lengtharray=arcpy.Array([point1.firstPoint,point2.firstPoint])
        LENGTH=arcpy.Polyline(lengtharray)

        cursorlength = arcpy.da.InsertCursor(cirque_length, ["SHAPE@", "ID_cirque"])
        cursorlength.insertRow([LENGTH, row[4]])
        del cursorlength

        #populate field in cirque shapefile
        row[0]=int(LENGTH.length)

        #width calculation
        crosspoint=LENGTH.centroid #center of the length line
        angle2=-(1/slope) #inverse of the tangent
        angleok=math.atan(angle2)#this is the perpendicular angle to the line
        pointxlength=(math.cos(angleok)*(LENGTH.length))*10 # direction and length of x coordinate for cross line from the center at one side
        #we assume that no cirque's length is ten times larger than width
        pointXpos=crosspoint.X+pointxlength # position of x at one side of the cross line
        pointylength=(math.sin(angleok)*(LENGTH.length))*10  # direction and length of y coordinate for cross line from the center at one side
        pointYpos=crosspoint.Y+pointylength # position of y at one side of the cross line
        pointA=arcpy.Point(pointXpos,pointYpos) #first corner of cross line
        pointXneg=crosspoint.X-pointxlength # position of x at other side of the cross line
        pointYneg=crosspoint.Y-pointylength # position of y at other side of the cross line
        pointB=arcpy.Point(pointXneg, pointYneg) #second corner of cross line

        widtharray1=arcpy.Array([pointA,crosspoint]) #first cross line
        WIDTH1=arcpy.Geometry("Polyline",widtharray1, spatialref) #first cross line geometry
        widtharray2=arcpy.Array([pointB,crosspoint]) #second cross line
        WIDTH2=arcpy.Geometry("Polyline",widtharray2, spatialref) #second cross line geometry

        #turn boundary into geometry
        array=bound.getPart()
        boundary=arcpy.Geometry("Polyline",array,spatialref)

        point1=WIDTH1.intersect(boundary,1) # first intersection of boundary and cross line
        point2=WIDTH2.intersect(boundary,1) # second intersection of boundary and cross line

        #width line
        widtharray=arcpy.Array([point1.firstPoint,point2.firstPoint])
        WIDTH=arcpy.Polyline(widtharray)
        cursorwidth = arcpy.da.InsertCursor(cirque_width, ["SHAPE@", "ID_cirque"])
        cursorwidth.insertRow([WIDTH, row[4]])
        del cursorwidth

        #populate field in cirque shapefile
        row[1]=int(WIDTH.length)

        #calculate L/W
        row[2]=row[0]/row[1]

        #update rows content so that it is permanently stored
        cursor.updateRow(row)

        #else:
        #    arcpy.AddMessage("cirque %r does not have a threshold or has two or more thresholds (e.g. is part of two or more nested cirques)" %(row[4])) #this relates to the cases when a cirque might not have a threshold point inside

del row, cursor


arcpy.CopyFeatures_management(cirque_length, OutputLength)
arcpy.CopyFeatures_management(cirque_width, OutputWidth)
arcpy.Delete_management("in_memory") ### Empty the in_memory



##Step 4: derive 3d statistics and hypsometry 
arcpy.AddMessage("Step 4: derive 3d statistics and hypsometry")

midAltContur = arcpy.CreateFeatureclass_management("in_memory", "midAltContur", "POLYLINE", "","","",cirques_copy)

#FcID = arcpy.Describe(cirques_copy).OIDFieldName
FcID = "ID_cirque"

fields = ("SHAPE@", "Z_min","Z_max","H","Z_mean","A3D","Slope_mean", "Aspectmean", "Plan_closISE", "Z_mid", "A3D_A2D", "Hypsomax", "HI","Prof_clos", FcID, "Z_median", "Asp_east","Asp_north", "Slope_max", "Slope_min", "Slpgt33","Slplt20","Slp20to33", "Plan_closSPA", "A2D")
volumetable = arcpy.env.scratchFolder + "\\volumetable.txt"
contur = arcpy.env.scratchGDB + "\\contur"
#startline = arcpy.env.scratchGDB + "\\startline"
#endline = arcpy.env.scratchGDB + "\\endline"
##it is better to calcuate the slope and aspect for the whole DEM first to keep the values for the whole outline
DEM_slope = Slope(InputDEM)  
DEM_aspect = Aspect(InputDEM)

with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
    i = 0
    for row in cursor:
        cirqueDTM = ExtractByMask(InputDEM, row[0]) ##shape@
        #cirqueSLOPE = Slope(cirqueDTM)  ##it is better to calcuate the slope and aspect for the whole DEM first to keep the values for the whole outline
        #cirqueASPECT = Aspect(cirqueDTM)
        cirqueSLOPE = ExtractByMask(DEM_slope, row[0])  ##it is better to calcuate the slope and aspect for the whole DEM first to keep the values for the whole outline
        cirqueASPECT = ExtractByMask(DEM_aspect, row[0])
        cirqueASPECT_rad = (cirqueASPECT * math.pi / 180)
        cirqueASPECT_sin = Sin (cirqueASPECT_rad)
        cirqueASPECT_cos = Cos (cirqueASPECT_rad)
        #DTM_min = arcpy.GetRasterProperties_management(cirqueDTM, "MINIMUM")
        #row[1] = int(float(DTM_min.getOutput(0)))
        row[1] = int(cirqueDTM.minimum)
        #DTM_max = arcpy.GetRasterProperties_management(cirqueDTM, "MAXIMUM")
        #row[2] = int(float(DTM_max.getOutput(0)))
        row[2] = int(cirqueDTM.maximum)
        row[3] = (row[2]-row[1])
        row[9] = (row[2]+row[1])/2 ##Middle elevation
        #DTM_mean = arcpy.GetRasterProperties_management(cirqueDTM, "MEAN")
        #row[4] = DTM_mean.getOutput(0)
        row[4] = cirqueDTM.mean
        
        #SLOPE_mean = arcpy.GetRasterProperties_management(cirqueSLOPE, "MEAN")
        #row[6] = SLOPE_mean.getOutput(0)
        row[6] = cirqueSLOPE.mean
        row[18] = cirqueSLOPE.maximum
        row[19] = cirqueSLOPE.minimum
        
        #SLOPE_max = arcpy.GetRasterProperties_management(cirqueSLOPE, "MAXIMUM")
        #SLOPE_min = arcpy.GetRasterProperties_management(cirqueSLOPE, "MINIMUM")
        #row[13] = float(SLOPE_max.getOutput(0)) - float(SLOPE_min.getOutput(0))
        row[13] = float(cirqueSLOPE.maximum) - float(cirqueSLOPE.minimum)
        
        #ASPECT_sin_mean = arcpy.GetRasterProperties_management(cirqueASPECT_sin, "MEAN")
        #ASPECT_sin_mean_value = float(ASPECT_sin_mean.getOutput(0))
        ASPECT_sin_mean_value = cirqueASPECT_sin.mean
        
        #ASPECT_cos_mean = arcpy.GetRasterProperties_management(cirqueASPECT_cos, "MEAN")
        #ASPECT_cos_mean_value = float(ASPECT_cos_mean.getOutput(0))
        ASPECT_cos_mean_value = cirqueASPECT_cos.mean
        
        radians_aspect = math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value)
        if  ASPECT_sin_mean_value  > 0 and ASPECT_cos_mean_value > 0 :
            degree_aspect = radians_aspect * 180/math.pi
            row[7] = degree_aspect
            sin_radians_aspect = math.sin(radians_aspect)
            cos_radians_aspect = math.cos(radians_aspect)
            #row [7] = float((math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value))*180/math.pi)
        elif ASPECT_cos_mean_value < 0 :
            degree_aspect = (radians_aspect + math.pi) * 180/math.pi
            row[7] = degree_aspect
            sin_radians_aspect = math.sin(radians_aspect + math.pi)
            cos_radians_aspect = math.cos(radians_aspect + math.pi)
            
            #row [7] = float(((math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value))*180/math.pi)+180)
        else:
            degree_aspect = (radians_aspect + 2* math.pi) * 180/math.pi
            row[7] = degree_aspect
            sin_radians_aspect = math.sin(radians_aspect + 2 * math.pi)
            cos_radians_aspect = math.cos(radians_aspect + 2 * math.pi)
            
            #row [7] = float(((math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value))*180/math.pi)+360)
        #arcpy.AddMessage("Aspect_sin_mean new: " + str(sin_radians_aspect))
        #arcpy.AddMessage("Aspect_cos_mean new: " + str(cos_radians_aspect))
        row[16] = sin_radians_aspect
        row[17] = cos_radians_aspect

        try:
            result = plan_closISE(cirqueDTM, contur) ####, startline, endline)
            planclosISE = result[0] ##planclos
            #row[3] = result[1] ##startAngle
            #row[4] = result[2] ##EndAngle
            #row[5] = result[3] ##AngType
            #row[6] = result[4] ##Flag
            row[8]= planclosISE
            arcpy.Append_management(contur, midAltContur, "NO_TEST")
        except:
            arcpy.AddMessage("There is an error to derive the plan_closISE!")
            pass            

        ##Plan_closSPA
        try:
            angle = plan_closSPA(cirqueDTM)
            row[23]= angle
        except:
            arcpy.AddMessage("There is an error to derive the plan_closSPA!")
            pass            
        
        #calculate 3D surface
        ##Step 1: Conduct Suface Volume analysis to generate the surface volume table, volumetable
        if arcpy.Exists(volumetable):
            arcpy.Delete_management(volumetable)

        min_elev = cirqueDTM.minimum
        #arcpy.AddMessage(min_elev)
        int_min_elev = int (min_elev - 1)
        #arcpy.SurfaceVolume_3d(cirqueDTM, volumetable, "ABOVE", "0")
        arcpy.SurfaceVolume_3d(cirqueDTM, volumetable, "ABOVE", str(int_min_elev))

        ##Step 2: Read the volume table for 3D area and 2Darea, and calculate the A3D/A2D ratio
        arr=arcpy.da.TableToNumPyArray(volumetable, ('AREA_2D', 'AREA_3D'))
        area_2D = float(arr[0][0])
        area_3D = float(arr[0][1])
        #Ratio3D2D = area_3D / area_2D
        if area_2D > 0:
            Ratio3D2D = area_3D / area_2D
        else: ## if area_2D is zero (very small polygons), set the ratio3d2d as 1
            Ratio3D2D = 1

        ##Step 3: Assign the values to the attibute fields
        row[5] = area_3D
        #row[9] = area_2D
        row[10] = Ratio3D2D
        ##Step 4: make sure to delete volumetable, so that it only has one record for the corresponding cirque outline
        arcpy.Delete_management(volumetable)

        ##Calcualte the Hypsometric max and Hypsometric intergal
        array = arcpy.RasterToNumPyArray(cirqueDTM,"","","",int_min_elev)
        EleArr = array[array > int_min_elev].astype(int) ##Get the elevations greater than zero
        Z_min = np.min(EleArr)
        Z_max = np.max(EleArr)
        Z_mean = np.mean(EleArr)
        Z_median = np.median(EleArr)
        row[15] = Z_median 

        Hi = (Z_mean - Z_min) / (Z_max - Z_min)
        row[12] = Hi
        
        vals,counts = np.unique(EleArr, return_counts=True)
        index = np.argmax(counts)
        hypo_max = vals[index]
        row[11] = hypo_max

        ##Derive the three added slope values: Slpgt33,Slplt20,Slp20to33
        array = arcpy.RasterToNumPyArray(cirqueSLOPE,"","","",-1)
        slpArr = array[array > -1].astype(float) ##Get all slope cells 
        total = len(slpArr)
        #slpgt33Arr = array[array > 33].astype(int) ##Get the slope greater than 33 degree
        slpgt33Arr = slpArr[slpArr > 33].astype(float) ##Get the slope greater than 33 degree
        slpgt33 = len(slpgt33Arr)
        slplt20Arr = slpArr[slpArr < 20].astype(float) ##Get the slope less than 20 degree
        slplt20 = len(slplt20Arr)
        row[20] = slpgt33 / total * 100 ##the percent of slope > 33 degree 
        row[21] = slplt20 / total * 100 ##the percent of slope < 20 degree 
        row[22] = (total - slpgt33 - slplt20) / total * 100        ##the percent of 20< slope < 33 degree

        #update cursor
        cursor.updateRow(row)
        i += 1

        arcpy.AddMessage("Finished cirque #" + str(i))

#delete cursor variables        
del row, cursor

##Derive, CS, L/H and W_H
fields = ("L", "W", "H", "CS", "L_H", "W_H")

with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
    for row in cursor:
        l = row[0]
        w = row[1]
        h = row[2]
        row[3] = pow(l*w*h, 1.0/3.0)
        row[4] = l / h
        row[5] = w / h

        cursor.updateRow(row)

#delete cursor variables        
del row, cursor


arcpy.CopyFeatures_management(midAltContur, OutputMidAltContour)
arcpy.Delete_management("in_memory") ### Empty the in_memory
arcpy.Delete_management(contur)
#arcpy.Delete_management(startline)
#arcpy.Delete_management(endline)


##Step 5: derive New cirque metrics info 
arcpy.AddMessage("Step 5: derive axis-related cirque metrics: Axprofileclos, Ax intergal, Ax Aspect, Amplitude, Axgrad, Length curve fitting parameters, and Width curve fitting parameters")

selLength = "in_memory\\selLength"
selWidth = "in_memory\\selWidth"

##Check the direction and flip the length from low to high elevations
arcpy.AddMessage("----derive length axis-related cirque metrics and curve fitting parameters")
Check_If_Flip_Line_Direction(OutputLength, InputDEM)
arcpy.InterpolateShape_3d(InputDEM, OutputLength, "in_memory\\length3D", 90) ##smoothed using 3 cell-size 
arcpy.InterpolateShape_3d(InputDEM, OutputWidth, "in_memory\\width3D", 90) ##smoothed using 3 cell-size 

FID_list = []
HLHI_list = []
HLAsp_list = []
P_clos_list = []
Amplitude_list = []
Axgrad_list = []

L_Exp_A_list = []
L_Exp_B_list = []
L_Exp_R2_list = []
L_Kcurv_C_list = []
L_Kcurv_R2_list = []
W_Quad_C_list = []
W_Quad_R2_list = []

with arcpy.da.SearchCursor("in_memory\\length3D", ["ID_Cirque","SHAPE@", "SHAPE@Length"]) as cursor:
    i = 0
    for row in cursor: ##Loop for each line
        PointX = []
        PointY = []
        LengthfromStart = []
        PointZ = []
        FID_list.append(row[0])
        lineLength = float(row[2])
        cumLength = 0
        for part in row[1]:
            pntCount = 0
            cumLength = 0
            segmentdis = 0
            for pnt in part:
                if pnt:
                    if pntCount > 0:
                        cumLength += Dist(startx, starty, pnt.X, pnt.Y) 
                    PointX.append(pnt.X)
                    PointY.append(pnt.Y)
                    PointZ.append(pnt.Z)
                    LengthfromStart.append(cumLength)

                    startx = pnt.X
                    starty = pnt.Y
                    pntCount += 1
        ##Calculate the HI value
        max_Z = max(PointZ)
        min_Z = min(PointZ)
        mean_Z = sum(PointZ) / len(PointZ)
        HI = (mean_Z - min_Z) / (max_Z - min_Z)+ 0.001 ##add 0.001 to avoid the divide of zero
        #arcpy.AddMessage("High-length HI: " + str(HI))
        HLHI_list.append(HI)
        Amplitude_list.append(max_Z - min_Z)
        gradient = 180.0/math.pi * math.atan((max_Z - min_Z)/max(LengthfromStart))
        #arcpy.AddMessage("Axial gradient: " + str(gradient))
        Axgrad_list.append(gradient)

        ##Calculate the HL-Aspect
        #dz  = PointZ[-1] - PointZ[0]
        dx  = PointX[0] - PointX[-1]
        dy  = PointY[0] - PointY[-1]
        #arcpy.AddMessage(str(dx))
        #arcpy.AddMessage(str(dy))

        aspect = 180.0/math.pi * math.atan2(dy, dx)
        #arcpy.AddMessage("Aspect is: " + str(aspect))
        if aspect < 90:
            adj_aspect = 90.0 - aspect
        else:
            adj_aspect = 360 + 90.0 - aspect

        #arcpy.AddMessage("Aspect is: " + str(adj_aspect))
        #row[5] = adj_aspect
        HLAsp_list.append(adj_aspect)

        ##Derive the exponential model fit for the longtitude profile 03/15/2023
        pointH = [y - min_Z for y in PointZ]
        #arcpy.AddMessage(pointH)
        #arcpy.AddMessage(LengthfromStart)

        HArr = np.array(pointH)
        LenArr = np.array(LengthfromStart)

        validHArr = HArr[HArr > 0]
        validLenArr = LenArr[HArr > 0]
        
        #expfit_results = exp_curve_fit(LengthfromStart[1:], pointH[1:])
        try:
            #polyfit_results = polyfit(validLenArr, [math.log(y) for y in pointH[1:]], 1)
            polyfit_results = polyfit(validLenArr, np.log(validHArr), 1)
            #a = polyfit[0]
            #arcpy.AddMessage(polyfit_results)
            b = polyfit_results['polynomial'][0]
            a = np.exp(polyfit_results['polynomial'][1])
            R2 = polyfit_results['determination']

            #arcpy.AddMessage("Exp_fit b is: " + str(b))       
            #arcpy.AddMessage("Exp_fit a is: " + str(a))       
            #arcpy.AddMessage("Exp_fit R2 is: " + str(R2))
        except:
            arcpy.AddMessage("There is an error!")
            b = -999
            a = -999
            R2 = -999

        L_Exp_A_list.append(a)
        L_Exp_B_list.append(b)
        L_Exp_R2_list.append(R2)
        
        ###Calculate the profile closure
        #arcpy.AddMessage(LengthfromStart)
        #arcpy.AddMessage(PointZ)
        startx = np.array(LengthfromStart[0:-1])
        endx = np.array(LengthfromStart[1:])
        startz = np.array(PointZ[0:-1])
        endz = np.array(PointZ[1:])
        dzdx = (endz - startz)/(endx - startx)

        slopes = 180/np.pi * np.arctan(dzdx)

        #arcpy.AddMessage(slopes)
        if len(slopes) > 3:
            min_slp = np.min(slopes[0:3])
            max_slp = np.max(slopes[-3:])
        else:
            min_slp = np.min(slopes)
            max_slp = np.max(slopes)

        p_close = max_slp - min_slp
        #arcpy.AddMessage("the profile closure is: " + str(p_close))
        P_clos_list.append(p_close)

        #K-curve-fit
        max_len = max(LengthfromStart)
        PointZ.reverse()
        LengthfromStart.reverse()
        normalH = np.array([(y - min_Z)/(max_Z - min_Z) for y in PointZ])
        normalLen = np.array([(max_len - y) /(max_len) for y in LengthfromStart])
        #arcpy.AddMessage(normalLen)
        #arcpy.AddMessage(normalH)
        
        fit_results = k_curve_fit(normalLen, normalH)
        c = fit_results[0]
        R2 = fit_results[1]

        L_Kcurv_C_list.append(c)
        L_Kcurv_R2_list.append(R2)
        
        i += 1
        arcpy.AddMessage("Finished cirque #" + str(i))

#delete cursor variables        
del row, cursor

##Do the width regression
arcpy.AddMessage("----derive quadratic fitting parameters along the cirque width axis")
                         
with arcpy.da.SearchCursor("in_memory\\width3D", ["ID_Cirque","SHAPE@"]) as cursor:
    i = 0
    for row in cursor: ##Loop for each line
        PointX = []
        PointY = []
        LengthfromStart = []
        PointZ = []
        #FID_list.append(row[0])
        #lineLength = float(row[2])
        cumLength = 0
        for part in row[1]:
            pntCount = 0
            cumLength = 0
            segmentdis = 0
            for pnt in part:
                if pnt:
                    if pntCount > 0:
                        cumLength += Dist(startx, starty, pnt.X, pnt.Y) 
                    PointX.append(pnt.X)
                    PointY.append(pnt.Y)
                    PointZ.append(pnt.Z)
                    LengthfromStart.append(cumLength)

                    startx = pnt.X
                    starty = pnt.Y
                    pntCount += 1

        ##Derive quadratic equation fit for the cross section profile along the width line 03/15/2023
        polyfit_results = polyfit(LengthfromStart,PointZ, 2)
        #a = polyfit[0]
        #arcpy.AddMessage(polyfit_results)
        c = polyfit_results['polynomial'][0]
        R2 = polyfit_results['determination']
        #arcpy.AddMessage("a is: " + str(a))       
        #arcpy.AddMessage("R2 is: " + str(R2))       
        W_Quad_C_list.append(c)
        W_Quad_R2_list.append(R2)         

        i += 1
        arcpy.AddMessage("Finished cirque #" + str(i))

#delete cursor variables        
del row, cursor


##Assign values to cirque polygons
#arcpy.AddMessage(FID_list)
#FcID = arcpy.Describe(cirques_copy).OIDFieldName
FcID = "ID_cirque"

fields = (FcID, "Axprofclos", "Axhli", "Axasp", "Axamp", "Axgrad", "L_Exp_A","L_Exp_B","L_Exp_R2","L_Kcurv_C","L_Kcurv_R2","W_Quad_C", "W_Quad_R2")

with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
    for row in cursor:
        try:
            fid = FID_list.index(row[0])
            #arcpy.AddMessage("fid is " + str(fid))
            row[1] = P_clos_list[fid]
            row[2] = HLHI_list[fid]
            row[3] = HLAsp_list[fid]
            row[4] = Amplitude_list[fid]
            row[5] = Axgrad_list[fid]

            row[6] = L_Exp_A_list[fid]
            row[7] = L_Exp_B_list[fid]
            row[8] = L_Exp_R2_list[fid]
            row[9] = L_Kcurv_C_list[fid]
            row[10] = L_Kcurv_R2_list[fid]
            row[11] = W_Quad_C_list[fid]
            row[12] = W_Quad_R2_list[fid]

            #update cursor
            cursor.updateRow(row)
        except:
            arcpy.AddMessage("There is an error in the calculation. Move to the next one")
            pass
           

#delete cursor variables        
del row, cursor


arcpy.Delete_management("in_memory") ### Empty the in_memory

##Step 6: derive New cirque metrics info 
arcpy.AddMessage("Step 6: derive catchment-related cirque metrics: Maxabalt and Pctabarea")

##Derive the upstream catchment for each cirque threshold
arcpy.AddMessage("----Deriving the upstream catchment for each cirque")

fillDEM =Fill(InputDEM)  ##Fill the sink first
fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
facc = FlowAccumulation(fdir) ##Flow accmulation

##Derive the streams and cross sections for cirque thresholds
Singlepoint = "in_memory\\Singlepoint"
Streams = "in_memory\\Streams"
FCselected = "in_memory\\FCselected"  
Singlestream = "in_memory\\Singlestream"
SingleWs = "in_memory\\SingleWs"
#singleBND = "in_memory\\singleBND"

#Get the minimum area of the cirque dataset
cirqueArr = arcpy.da.FeatureClassToNumPyArray(cirques_copy, "SHAPE@AREA")
areas = np.array([item[0] for item in cirqueArr])
min_cirque_area = np.min(areas)
#arcpy.AddMessage(str(min_cirque_area))

threshold_cells = int(min_cirque_area / (cellsize_int * cellsize_int)/4) ##divde 2 to make sure there is stream out from the cirque
#arcpy.AddMessage(str(threshold_cells))
threshold = min(50, threshold_cells)
#arcpy.AddMessage(str(threshold))

outGreaterThan = Con(facc > threshold, 1,0)  ##Determine the highest flowaccumuation part
# Process: Stream Link
outStreamLink = StreamLink(outGreaterThan, fdir)
  
# Process: Stream to Feature
StreamToFeature(outStreamLink, fdir, Streams, "NO_SIMPLIFY")


MaxFccTable = "in_memory\\MaxFccTable"
TmpStream = "in_memory\\TmpStream"

StreamToFeature(outStreamLink, fdir, TmpStream, "SIMPLIFY")
ZonalStatisticsAsTable(outStreamLink, "VALUE", facc, MaxFccTable, "DATA", "MAXIMUM")
# Process: Join Field
arcpy.JoinField_management(TmpStream, "grid_code", MaxFccTable, "Value", "MAX")  ##Join to get the flow accumulation value

##Convert cirque outline polygon to polylines
cirque_lines = "in_memory\\cirque_lines"
arcpy.PolygonToLine_management(cirques_copy, cirque_lines, "IGNORE_NEIGHBORS")

##Derive the watershed for cirque thresholds
TotalWS = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, "TotalWS","POLYGON","","","",spatialref)
arcpy.AddField_management(TotalWS, "ID_cirque", "INTEGER", 6)

countResult = arcpy.GetCount_management(cirque_lines)
count = int(countResult.getOutput(0))
#FcID = arcpy.Describe(cirque_lines).OIDFieldName
#FcID = arcpy.Describe(inputFCcopy).OIDFieldName
FcID = "ID_cirque"

linearr = arcpy.da.FeatureClassToNumPyArray(cirque_lines, FcID)
id_list = np.array([item[0] for item in linearr])

#arcpy.AddMessage(id_list)

for fid in id_list:
    #arcpy.AddMessage("Generating cirque "+str(ifeature + 1)+" of "+str(count))
    arcpy.AddMessage("Processing cirque #" + str(fid))
    query = FcID +" = "+str(fid)
    #arcpy.AddMessage(query)
    arcpy.Select_analysis(cirque_lines, FCselected, query)
    
    arcpy.Intersect_analysis([Streams, FCselected], Singlepoint, "#", "#", "POINT")

    pntcountResult = arcpy.GetCount_management(Singlepoint)
    pntcount = int(pntcountResult.getOutput(0))
    #arcpy.AddMessage("the number of point is:" + str(pntcount))
    if (pntcount == 0): ##if no intersect points, use the buffer to get the intersection points for another time
        tmpbuf = "in_memory\\tmpbuf"
        arcpy.Buffer_analysis(FCselected, tmpbuf, "5 Meters")
        arcpy.Intersect_analysis([Streams, tmpbuf], Singlepoint, "#", "#", "POINT")
        pntcountResult = arcpy.GetCount_management(Singlepoint)
        pntcount = int(pntcountResult.getOutput(0))
        #arcpy.AddMessage("the number of point is:" + str(pntcount))
        arcpy.Delete_management (tmpbuf)
        
    if (pntcount > 0):
        #Calculate Watershed
        outPour = SnapPourPoint(Singlepoint, facc, 0)
        outWs1 = Watershed(fdir, outPour)
        outWs = Con(outWs1 >= 0, 1)  
        arcpy.RasterToPolygon_conversion(outWs, SingleWs, "NO_SIMPLIFY", "VALUE")

        ##Append the cirque polygon to the singleWS to prevent very small ws created: Added by Yingkui Li 10/26/2023
        arcpy.Select_analysis(cirques_copy, "in_memory\\Selcirque", query)
        arcpy.Append_management("in_memory\\Selcirque", SingleWs, "NO_TEST")
         
        arcpy.Dissolve_management(SingleWs, "in_memory\\dissolve_SingleWs")
        arcpy.AddField_management("in_memory\\dissolve_SingleWs", "ID_cirque", "INTEGER", 6)
        with arcpy.da.UpdateCursor("in_memory\\dissolve_SingleWs", "ID_cirque") as cursor:
            for row in cursor:
               row[0]=fid
               cursor.updateRow(row)
        del row, cursor
        arcpy.Append_management("in_memory\\dissolve_SingleWs", TotalWS, "NO_TEST")        

    

#arcpy.CopyFeatures_management(cirque_thresholds, OutputThresholds)

arcpy.Delete_management("in_memory") ### Empty the in_memory

#arcpy.CopyFeatures_management(TotalWS, "d:\\temp\\TotalWS.shp")

##Derive the upstream catchment for each cirque threshold
arcpy.AddMessage("----Deriving Maxabalt and Pctabarea")

#FcID = arcpy.Describe(cirques_copy).OIDFieldName
FcID = "ID_cirque"

fields = ("SHAPE@", "SHAPE@AREA","Z_max","Z_mean", "Maxabalt", "Pctabarea", FcID)

selWs = "in_memory\\selWs"

with arcpy.da.UpdateCursor(cirques_copy, fields) as cursor:
    i = 0
    for row in cursor:
        query = "ID_cirque = "+ str(row[6])
        arcpy.Select_analysis(TotalWS, selWs, query)
        
        arr = arcpy.da.FeatureClassToNumPyArray(selWs, "SHAPE@AREA")
        areas = np.array([item[0] for item in arr])
        #arcpy.AddMessage(str(len(areas)))
        area = max(areas) ##only choose the maximum area
        #arcpy.AddMessage(str(area))
        WsDTM = ExtractByMask(InputDEM, selWs) ##shape@
        #WsDTM_max = arcpy.GetRasterProperties_management(WsDTM, "MAXIMUM")
        #maxabalt = int(float(WsDTM_max.getOutput(0)))
        maxabalt = int(float(WsDTM.maximum))
        if maxabalt < row[2]: ##in case the maxabalt is smaller than the Z_max from the cirque outline
            row[4] = row[2]
            row[5] = 100
        else:
            mean_ele = row[3] 
            row[4] = maxabalt
            ##Need to figure out the area above the cirques
            ##Erase the watershed using cirque polygon
            arcpy.Erase_analysis(selWs, row[0], "in_memory\\ws_leftover")
            wsarr = arcpy.da.FeatureClassToNumPyArray("in_memory\\ws_leftover", "OID@")
            ab_area = 0
            if len(wsarr) > 0:
                ##multiple parts to single parts
                arcpy.MultipartToSinglepart_management("in_memory\\ws_leftover", "in_memory\\ws_leftover_singleParts")
                ##delete small spurious polygons *1 cell sizes
                with arcpy.da.UpdateCursor("in_memory\\ws_leftover_singleParts", 'SHAPE@AREA') as cursor2:  ##outline_cp is the outmost simplified outlines
                    for row2 in cursor2:
                        if row2[0] < 1000: ##about 1 cell-sizes pixcels
                            cursor2.deleteRow()
                del row2, cursor2
                wsarr2 = arcpy.da.FeatureClassToNumPyArray("in_memory\\ws_leftover_singleParts", "OID@")
                if len(wsarr2) > 0:
                    ##get the elevation info for the polygons
                    Zonaltable = arcpy.env.scratchGDB + "\\Zonaltable"
                    outZSaT = ZonalStatisticsAsTable("in_memory\\ws_leftover_singleParts", "ObjectID", WsDTM, Zonaltable, "NODATA", "MEAN")
                    arcpy.JoinField_management("in_memory\\ws_leftover_singleParts", "ObjectID", Zonaltable, "ObjectID_1", ["MEAN"])
                    arcpy.Delete_management(Zonaltable)
                    ##Get the total areas with mean > Z_mean 
                    with arcpy.da.SearchCursor("in_memory\\ws_leftover_singleParts", ['SHAPE@AREA', 'MEAN']) as cursor3:  ##outline_cp is the outmost simplified outlines
                        for row3 in cursor3:
                            #arcpy.AddMessage("Mean is: " + str(row3[1]))
                            try:
                                if row3[1] > mean_ele: 
                                    ab_area += row3[0]    
                            except:
                                #arcpy.AddMessage("Mean is: " + str(row3[1]))
                                pass
                    del row3, cursor3
            #arcpy.AddMessage(str(ab_area))
            row[5] = row[1]/(row[1]+ab_area)*100

        #update cursor
        cursor.updateRow(row)

        i += 1
        arcpy.AddMessage("Finished cirque #" + str(i))

#delete cursor variables        
del row, cursor

arcpy.CopyFeatures_management(cirques_copy, outputCirques)

arcpy.Delete_management(cirques_copy) 
arcpy.Delete_management(TotalWS) 

arcpy.Delete_management("in_memory") ### Empty the in_memory












