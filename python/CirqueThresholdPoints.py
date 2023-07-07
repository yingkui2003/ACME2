#-------------------------------------------------------------------------------
# Name:        CirqueThresholdPoints
# Purpose:     Derive Cirque threshold points based on the cirque outlines and DEM
# Author:      Yingkui Li
# Created:     07/19/2022 - 08/23/2022
# Copyright:   (c) Yingkui Li 2022
#-------------------------------------------------------------------------------

#import locale
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
#import itertools
import math
import os, sys
import matplotlib.pyplot as plt
from scipy import optimize

arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "0.01 Meters"
#locale.setlocale(locale.LC_ALL,"")#sets local settings to decimals

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


def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
  

#---------------------------------------------------------------------------------------------------------------
# This function calculates the distance between two points
#--------------------------------------------------------------------------------------------------------------- 
def Dist(x1,y1,x2,y2):
    return math.sqrt(math.pow(math.fabs(x1-x2),2)+math.pow(math.fabs(y1-y2),2))

def connect_major_line_sections(line):
    lineArray = arcpy.da.FeatureClassToNumPyArray(line,"SHAPE@LENGTH")
    lengths = np.array([item[0] for item in lineArray])
    if len(lengths) == 1:
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
        #arcpy.AddMessage("GAP distance: " + str(distArr[index]))
        ##Figure out the start and end location of the points
        start = int(index / pntCount)
        end = index % pntCount
        ##get the length of the start line
        start_line_length = lengths[start]
        #arcpy.AddMessage("start_line_length: " + str(start_line_length))
        end_line_length = lengths[end]
        #arcpy.AddMessage("start_line_length: " + str(end_line_length))
        
        ##only connect the line when the distance between the start or end points
        ##of two line section is less than the minimum length of the two lines
        if distArr[index] < min(start_line_length, end_line_length):
            #arcpy.AddMessage("Connect two lines")
            array = arcpy.Array([arcpy.Point(startpntX[start],startpntY[start]),arcpy.Point(endpntX[end], endpntY[end])])
            polyline = arcpy.Polyline(array)
            new_line_cursor.insertRow([polyline])

    del new_line_cursor

    arcpy.Append_management(connectionline, line, "NO_TEST")
    arcpy.Dissolve_management(line, "in_memory\\line_dissolve", "", "", "SINGLE_PART")
    ##check to delete unconnected lines
    lineArray = arcpy.da.FeatureClassToNumPyArray("in_memory\\line_dissolve",["SHAPE@LENGTH"])
    if len(lineArray) > 1:
        #arcpy.AddMessage("Delete unconnected small line sections")
        lengths = np.array([item[0] for item in lineArray])
        #arcpy.AddMessage(lengths)
        max_length = max(lengths)
        with arcpy.da.UpdateCursor("in_memory\\line_dissolve",("SHAPE@LENGTH")) as cursor:
            for row in cursor:
                if row[0] < max_length:
                    cursor.deleteRow()
        #delete cursors
        del row, cursor
    
    arcpy.CopyFeatures_management("in_memory\\line_dissolve", line)

    arcpy.Delete_management ("in_memory\\line_dissolve")
    arcpy.Delete_management (connectionline)
    
    return line

def CirqueThresholds (cirque_polys, dem): ##, overlap = False):
    #get the count of autopolys
    #make copies of the data to in-memory data
    cirquepolys = "in_memory\\cirquepolys"

    cirquePoints = "in_memory\\cirquePoints"

    arcpy.CopyFeatures_management(cirque_polys, cirquepolys)

    ##Check if the cirques have overlaps
    b_Overlap = False
    outlines = arcpy.CopyFeatures_management(cirquepolys, arcpy.Geometry())
    for i in range(len(outlines)):
        #arcpy.AddMessage(str(i))
        for j in range(len(outlines)):
            if j > i:
                b_disjoint = outlines[i].disjoint(outlines[j])
                if b_disjoint == False:
                    #arcpy.AddMessage("Joint: " + str(i) + " and " + str(j))
                    #ids.append(i)
                    b_Overlap = True
                    break
        if b_Overlap:
            break

    fillDEM =Fill(dem)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
    facc = FlowAccumulation(fdir) ##Flow accmulation
    
    #countResult = arcpy.GetCount_management(cirquepolys)
    #count = int(countResult.getOutput(0))
    FcID = arcpy.Describe(cirquepolys).OIDFieldName


    cirquePoints = arcpy.CreateFeatureclass_management("in_memory", "cirquePoints","POINT","","","",cirque_polys)

    if b_Overlap == False:
        ##if without cirque overlaps, can do the whole analysis without looping for each cirque polygon
        ##revised and tested by Yingkui on 7/15/2022 and sucessful!!! but some times one polygon have multiple points
        outZonalStatistics = ZonalStatistics(cirquepolys, FcID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part

        arcpy.RasterToPoint_conversion(outPour, cirquePoints, "VALUE")

        ##Check if one cirque produces multiple points by the same grid_code and mark the points so that the user can choose the correct one later
        ##Get the existing fields in cirquepolys
        exist_no_use_fields = [f.name for f in arcpy.ListFields(cirquepolys)] #List of current field names in outline layer

        arcpy.AddField_management(cirquepolys, 'ID_CirqueID', 'Long', 6) ##OFID is not FID, but it is related to FID
        arcpy.CalculateField_management(cirquepolys,"ID_Cirque",str("!"+str(arcpy.Describe(cirquepolys).OIDFieldName)+"!"),"PYTHON_9.3")

        joinedpoints = "in_memory\\joinedpoints"
        arcpy.SpatialJoin_analysis(cirquePoints, cirquepolys, joinedpoints, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
        
        arr = arcpy.da.FeatureClassToNumPyArray(joinedpoints, 'ID_Cirque')
        Cirque_IDs = np.array([item[0] for item in arr])
        unique_CirqueIDs, counts = np.unique(Cirque_IDs, return_counts=True)

        multipeIDs = unique_CirqueIDs[counts > 1] ##only find the duplicated codes
        #arcpy.AddMessage(multipeIDs)
        if len(multipeIDs) > 0: ##check if there is multiple points for some cirques
            arcpy.AddMessage("Some cirques include multiple points! These points are flagged as 1 in the attrbute table")
            arcpy.AddField_management(joinedpoints, 'Flag', 'Long', 6) ##OFID is not FID, but it is related to FID
            with arcpy.da.UpdateCursor(joinedpoints,("ID_Cirque","Flag")) as cursor:
                for row in cursor:
                    CID = row[0]
                    row[1] = 0
                    if CID in multipeIDs:
                       row[1] = 1
                    cursor.updateRow(row)
            #delete cursors
            del row, cursor

        ##Delete no_use fields and the fields created by the spatial join
        #arcpy.AddMessage(exist_no_use_fields[2:])
        if len(exist_no_use_fields[2:]) > 0: 
            arcpy.DeleteField_management(joinedpoints, exist_no_use_fields[2:])
        arcpy.DeleteField_management(joinedpoints, ["Join_Count", "TARGET_FID", "grid_code"])

        return joinedpoints

    else:
        sel_poly = "in_memory\\sel_poly"
        Singlepoint = "in_memory\\Singlepoint"
        ##Need to loop for each cirque polygon to make sure that each cirque polygon has a threshold points
        countResult = arcpy.GetCount_management(cirquepolys)
        count = int(countResult.getOutput(0))
        FcID = arcpy.Describe(cirquepolys).OIDFieldName
        polyarray = arcpy.da.FeatureClassToNumPyArray(cirquepolys,FcID)
        FcID_Arr = np.array([item[0] for item in polyarray])
        
        for fid in FcID_Arr:
            arcpy.AddMessage("Processing cirque #" + str(fid))
            query = FcID +" = "+ str(fid)  ##THe FID starts from 1 for the in_memory data; Which is different from shape (FID from 0) and feature class in the geodatabase (objectID from 1)
            #arcpy.AddMessage(query)
            arcpy.Select_analysis(cirquepolys, sel_poly, query)

            PolyID = arcpy.Describe(sel_poly).OIDFieldName
            
            outZonalStatistics = ZonalStatistics(sel_poly, PolyID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
            outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
            
            arcpy.RasterToPoint_conversion(outPour, Singlepoint, "VALUE")
            arcpy.AddField_management(Singlepoint, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID
            with arcpy.da.UpdateCursor(Singlepoint,("ID_cirque")) as cursor:
                for row in cursor:
                    row[0] = fid
                    cursor.updateRow(row)
            #delete cursors
            del row, cursor
            
            ##Make sure only one point produced, otherwise, add a new field with flag
            pntcountResult = arcpy.GetCount_management(Singlepoint)
            pntcount = int(pntcountResult.getOutput(0))
            if pntcount > 1:
                arcpy.AddMessage("Some cirques include multiple points! These points are flagged as 1 in the attrbute table")
                arcpy.AddField_management(Singlepoint, 'Flag', 'Long', 6) ##OFID is not FID, but it is related to FID
                with arcpy.da.UpdateCursor(Singlepoint,("Flag")) as cursor:
                    for row in cursor:
                        row[0] = 1
                        cursor.updateRow(row)
                #delete cursors
                del row, cursor
            ##Append the single point to 
            arcpy.Append_management(Singlepoint, cirquePoints, "NO_TEST") 

        ##delete the temp dataset
        arcpy.Delete_management (sel_poly)
        arcpy.Delete_management (Singlepoint)

        return cirquePoints


def CirqueThresholds_midpoints(InputCirques, InputDEM, percentile):
    ##This method determines the cirque threshold mid-points by the lower 10th percentile of the elevation of the cirque outlines
    sel_poly = "in_memory\\sel_poly"
    cirque_line = "in_memory\\cirque_line"
    tmpbuf = "in_memory\\tmpbuf"
    points = "in_memory\\points"
    points_with_Z = "in_memory\\points_with_Z"
    sel_points = "in_memory\\sel_points"
    tmpline = "in_memory\\tmpline"
    tmpline3D = "in_memory\\tmpline3D"
    midpoint = "in_memory\\midpoint"

    ##Obtain the cellsize of the DEM
    dem = arcpy.Raster(InputDEM)
    cellsize = dem.meanCellWidth
    #arcpy.AddMessage("Cell size: " + str(cellsize))

    cirquePoints = arcpy.CreateFeatureclass_management("in_memory", "cirquePoints","POINT","","","",InputCirques)
    arcpy.AddField_management(cirquePoints, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID

    cirqueThresholds = arcpy.CreateFeatureclass_management("in_memory", "cirqueThresholds","POLYLINE","","","",InputCirques)
    arcpy.AddField_management(cirqueThresholds, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID
    

    #FcID = "ORIG_FID"
    FcID = arcpy.Describe(InputCirques).OIDFieldName
    linearray = arcpy.da.FeatureClassToNumPyArray(InputCirques,FcID)
    FcID_Arr = np.array([item[0] for item in linearray])
    i = 0    
    for fid in FcID_Arr:
        arcpy.AddMessage("Processing cirque #" + str(fid))
        query = FcID +" = "+ str(fid)  ##THe FID starts from 1 for the in_memory data; Which is different from shape (FID from 0) and feature class in the geodatabase (objectID from 1)
        #arcpy.AddMessage(query)
        arcpy.Select_analysis(InputCirques, sel_poly, query)

        arcpy.PolygonToLine_management(sel_poly, cirque_line, "IGNORE_NEIGHBORS")    
        arcpy.FeatureVerticesToPoints_management(cirque_line, points, "ALL")

        arcpy.Buffer_analysis(sel_poly, tmpbuf, "100 Meters")
        cirqueDTM = ExtractByMask(dem, tmpbuf)
        
        ExtractValuesToPoints(points, dem, points_with_Z, "INTERPOLATE", "VALUE_ONLY")
        pntArray = arcpy.da.FeatureClassToNumPyArray(points_with_Z,["RASTERVALU"])
        Elevs = np.array([item[0] for item in pntArray])
        #arcpy.AddMessage(Elevs)

        cutoff_elev = np.percentile(Elevs, percentile) ##get the 10th percentile elevation
        #arcpy.AddMessage("Cutoff Elevation from outline: " + str(cutoff_elev))

        ##Method 1: use the raster functions to get the cirque outlines within the lower elevations
        outCon = Con(cirqueDTM < cutoff_elev, 1)
        OutBndCln = BoundaryClean(outCon)
        arcpy.conversion.RasterToPolygon(OutBndCln, "in_memory\\OutBndCln_poly")
        arcpy.analysis.Clip(cirque_line, "in_memory\\OutBndCln_poly", tmpline)
        arcpy.MultipartToSinglepart_management(tmpline, "in_memory\\tmpline_singlePart")
        arcpy.Dissolve_management("in_memory\\tmpline_singlePart", tmpline, "", "", "SINGLE_PART")
        
        ##test if there are more line sections
        lineArray = arcpy.da.FeatureClassToNumPyArray(tmpline,["SHAPE@LENGTH"])
        if len(lineArray) > 1:
            #arcpy.AddMessage("multiple lines created, need to connect the lines")
            ##only connected major lines
            connect_major_line_sections(tmpline)
        
        ##Check the slopes along the line
        ##use filled DEM and 3*cellsize as spacing; save the 3d feature as one output: out3DProfiles
        arcpy.InterpolateShape_3d(dem, tmpline, tmpline3D, cellsize*3)
        arcpy.FeatureVerticesToPoints_management(tmpline3D, points, "ALL")
        pntArray = arcpy.da.FeatureClassToNumPyArray(points,["SHAPE@X", "SHAPE@Y", "SHAPE@Z"])
        pntX = np.array([item[0] for item in pntArray])
        pntY = np.array([item[1] for item in pntArray])
        pntZ = np.array([item[2] for item in pntArray])
        #arcpy.AddMessage(pntZ)

        maxZ = max(pntZ)
        minZ = min(pntZ)
        Z_range = maxZ - minZ
        Z_interval = 5 ##5 m interval to divide Z
        numZ = int(Z_range / Z_interval) 
        widthList = []
        ZList = []
        for i in range(numZ):
            Z = maxZ - i * Z_interval
            selX = pntX[pntZ < Z]
            selY = pntY[pntZ < Z]
            width = Dist(selX[0], selY[0], selX[-1], selY[-1])
            widthList.append(width)
            ZList.append(Z)
        
        #arcpy.AddMessage(ZList)
        #arcpy.AddMessage(widthList)

        x = np.array(ZList)
        y = np.array(widthList)
        try:
            p , e = optimize.curve_fit(piecewise_linear, x, y)
            #arcpy.AddMessage(p)
            #arcpy.AddMessage(e[0][0])
            '''
            ##create the x-y plot with the regression lines
            xd = np.linspace(minZ, maxZ, 100)
            #plt.plot(x, y, "o")
            #plt.plot(xd, piecewise_linear(xd, *p))
            #plt.show()
            '''
            if str(e[0][0]) != "inf":
                ##Use the elevation turning point to determine the threshold
                turn_Z = p[0]
                ##Find the point X and Y with the elevation < turn_Z+1
                selX = pntX[pntZ < turn_Z+1]
                selY = pntY[pntZ < turn_Z+1]
                #selZ = pntZ[pntZ < turn_Z+1]
                ##turn the XY points into a line feature
                #XYArr = np.column_stack([selX, selY, selZ])
                #arcpy.AddMessage(feature_info)
                #field = "Elev"
                startX = selX[0:-1]
                endX = selX[1:]
                startY = selY[0:-1]
                endY = selY[1:]

                new_line = arcpy.CreateFeatureclass_management("in_memory", "new_line","POLYLINE", "","","", tmpline)
                new_line_cursor = arcpy.da.InsertCursor(new_line, ('SHAPE@'))

                for i in range(len(startX)):
                    array = arcpy.Array([arcpy.Point(startX[i],startY[i]),arcpy.Point(endX[i], endY[i])])
                    polyline = arcpy.Polyline(array)
                    new_line_cursor.insertRow([polyline])

                del new_line_cursor

                arcpy.Dissolve_management(new_line, tmpline, "", "", "SINGLE_PART")
                arcpy.Delete_management(new_line)    
        except:
            #arcpy.AddMessage("Error")
            pass
        
        ##Get the center points of the tmpline
        arcpy.FeatureVerticesToPoints_management(tmpline, midpoint, "MID")

        arcpy.AddField_management(midpoint, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID
        with arcpy.da.UpdateCursor(midpoint,("ID_cirque")) as cursor:
            for row in cursor:
                row[0] = fid
                cursor.updateRow(row)
        #delete cursors
        del row, cursor

        arcpy.AddField_management(tmpline, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID
        with arcpy.da.UpdateCursor(tmpline,("ID_cirque")) as cursor:
            for row in cursor:
                row[0] = fid
                cursor.updateRow(row)
        #delete cursors
        del row, cursor

        
        ##Append the single point to 
        arcpy.Append_management(midpoint, cirquePoints, "NO_TEST") 
        arcpy.Append_management(tmpline, cirqueThresholds, "NO_TEST") 
        i += 1
        
    ##delete the temp dataset
    arcpy.Delete_management (sel_poly)
    arcpy.Delete_management (cirque_line)
    arcpy.Delete_management (tmpbuf)
    arcpy.Delete_management (points)
    arcpy.Delete_management (points_with_Z)
    arcpy.Delete_management (sel_points)
    arcpy.Delete_management (tmpline)
    arcpy.Delete_management (midpoint)
    arcpy.Delete_management ("in_memory\\OutBndCln_poly")
    arcpy.Delete_management ("in_memory\\OutBndCln_poly")
    arcpy.Delete_management ("in_memory\\tmpline_singlePart")
    
    return cirquePoints, cirqueThresholds


##--main program
##Main program
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)
method = arcpy.GetParameterAsText(2)
OutputCirques = arcpy.GetParameterAsText(3)

arcpy.Delete_management("in_memory") ### Empty the in_memory

#arcpy.AddMessage(method)

if method == "Mainstream exit":
    result_threshold_midpoints = CirqueThresholds (InputCirques, InputDEM)
else:
    result_threshold_midpoints, cirqueThresholds = CirqueThresholds_midpoints(InputCirques, InputDEM, 10)
    #arcpy.CopyFeatures_management(cirqueThresholds, "d:\\temp\\cirquethresholds.shp")
    
arcpy.CopyFeatures_management(result_threshold_midpoints, OutputCirques)

arcpy.Delete_management("in_memory") ### Empty the in_memory


