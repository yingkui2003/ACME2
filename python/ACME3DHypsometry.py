#-------------------------------------------------------------------------------
# Name:        Cirque 3D statistics and Hypsometry
# Purpose:     Derive the cirque metrics related to elevation, slope, and hypsometry
#
# Orignial Author:      Ramon Pellitero and Matteo Spagnolo
#
# Created:     30/04/2016
# Copyright:   (c) Ramon Pellitero and Matteo Spagnolo 2016
#
# Revised by Yingkui Li during 07/19/2022 - 8/24/2022
# The revised version combined the orginal 3D statistics and Hypsometry. 
# First, the A3D/A2D ratio is added because it may not be derived correctly using the Area3D
# derived from surface volume analysis and Area2D derived from the cirque polygon.
# Surface volume analysis is based on the TIN, and in some cases, the 3D area of the TIN is
# smaller than the Area2D derived from the polygon or raster. In the revised codes, the A3D/A2D
# ratio is directly derived using the 3D and 2D areas from the ArcGIS surface volume analysis.
# 
# Second, the original codes derive the hypsometric max and hypsometric intergal based on the contour
# line interval apporach, which is very time consuming. The revised codes derive the hypsometric metrics
# based on the mode of elevation values and the Elevation-relief ratio. The revised codes is much faster
# in determining the hypsometric max and hypsometric intergal values.
#
# The revised codes also improve the logic of the codes and fixed the memory allocation issue
# when processing a large amount of cirques.
#
# The revsied codes can work on both ArcMap and ArcGIS Pro.
#
# Created:     07/19/2022 - 08/23/2022
# Copyright:   (c) Yingkui Li 2022
#-------------------------------------------------------------------------------

from __future__ import division
import arcpy
from arcpy import env
from arcpy.sa import *
import math
import time
import numpy as np

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

#------------------------------------------------------------------------------------------------------------
# This function calculates the plan closure for a cirque. The codes are modified from ACME code
#------------------------------------------------------------------------------------------------------------
def plan_clos(cirqueDEM):
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

##Main program
# Script arguments
cirque = arcpy.GetParameterAsText(0)
dtm = arcpy.GetParameterAsText(1)

#environments
spatialref=arcpy.Describe(cirque).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

#a "list" where the name of fields from the attributed table are copied in
Fieldlist=[]
ListFields=arcpy.ListFields(cirque)
for x in ListFields:
    Fieldlist.append(x.baseName)

#if fields with the same name are already in the attribute table, no need to create them again (which will crash the tool)
fields = ("SHAPE@", "Z_min","Z_max","Z_range","Z_mean","area_3D","Slope_mean", "Aspectmean", "Plan_clos", "TIN_A2D", "A3d_A2D", "hypsomax", "HI")
#new_fields = ("Z_min","Z_max","Z_range","Z_mean","area_3D","Slope_mean", "Aspectmean", "Plan_clos", "TIN_A2D", "A3d_A2D", "hypsomax", "HI")
#arcpy.AddMessage(fields[1:])
for field in fields[1:]:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(cirque, field, "FLOAT",6 ,2)

volumetable = tempRs = arcpy.env.scratchFolder + "\\volumetable.txt"

with arcpy.da.UpdateCursor(cirque, fields) as cursor:
    i = 0
    for row in cursor:
        cirqueDTM = ExtractByMask(dtm, row[0]) ##shape@
        cirqueSLOPE = Slope(cirqueDTM)
        cirqueASPECT = Aspect(cirqueDTM)
        cirqueASPECT_rad = (cirqueASPECT * math.pi / 180)
        cirqueASPECT_sin = Sin (cirqueASPECT_rad)
        cirqueASPECT_cos = Cos (cirqueASPECT_rad)
        ##Change the following codes to avoid the use of comma as a decimal separator
        #DTM_min = arcpy.GetRasterProperties_management(cirqueDTM, "MINIMUM")
        #row[1] = int(float(DTM_min.getOutput(0)))
        row[1] = int(cirqueDTM.minimum)
        #DTM_max = arcpy.GetRasterProperties_management(cirqueDTM, "MAXIMUM")
        #row[2] = int(float(DTM_max.getOutput(0)))
        row[2] = int(cirqueDTM.maximum)
        row[3] = (row[2]-row[1])
        #DTM_mean = arcpy.GetRasterProperties_management(cirqueDTM, "MEAN")
        #row[4] = DTM_mean.getOutput(0)
        row[4] = cirqueDTM.mean

        #SLOPE_mean = arcpy.GetRasterProperties_management(cirqueSLOPE, "MEAN")
        #row[6] = SLOPE_mean.getOutput(0)
        row[6] = cirqueSLOPE.mean

        #ASPECT_sin_mean = arcpy.GetRasterProperties_management(cirqueASPECT_sin, "MEAN")
        #ASPECT_sin_mean_value = float(ASPECT_sin_mean.getOutput(0))
        ASPECT_sin_mean_value = cirqueASPECT_sin.mean
        
        #ASPECT_cos_mean = arcpy.GetRasterProperties_management(cirqueASPECT_cos, "MEAN")
        #ASPECT_cos_mean_value = float(ASPECT_cos_mean.getOutput(0))
        ASPECT_cos_mean_value = cirqueASPECT_cos.mean

        if  ASPECT_sin_mean_value  > 0 and ASPECT_cos_mean_value > 0 :
            row [7] = float((math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value))*180/math.pi)
        elif ASPECT_cos_mean_value < 0 :
            row [7] = float(((math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value))*180/math.pi)+180)
        else :
            row [7] = float(((math.atan (ASPECT_sin_mean_value/ASPECT_cos_mean_value))*180/math.pi)+360)

        #Get the plan_closure
        planclos = plan_clos(cirqueDTM)
        row[8]= planclos

        #calculate 3D surface
        ##The following codes are revised by Yingkui Li (7/19/2022) to derive the 3D area and A3D/A2D ratio based on ArcGIS surface volume analysis.
        ##The A3D/A2D ratio cannot be calcuated using the Area3D derived from surface volume analysis and Area2D derived from the cirque polygon
        ##because Surface volume analysis is based on the TIN, and the 3D area of the TIN is smaller than the Area2D derived from the polygon or raster
        ##Therefore, the A3D/A2D ratio is directly derived using the 3D and 2D areas from the ArcGIS surface volume analysis
        ##Step 1: Conduct Suface Volume analysis to generate the surface volume table, volumetable
        arcpy.SurfaceVolume_3d(cirqueDTM, volumetable, "ABOVE", "0")
        ##Step 2: Read the volume table for 3D area and 2Darea, and calculate the A3D/A2D ratio
        arr=arcpy.da.TableToNumPyArray(volumetable, ('AREA_2D', 'AREA_3D'))
        area_2D = float(arr[0][0])
        area_3D = float(arr[0][1])
        Ratio3D2D = area_3D / area_2D
        ##Step 3: Assign the values to the attibute fields
        row[5] = area_3D
        row[9] = area_2D
        row[10] = Ratio3D2D
        ##Step 4: make sure to delete volumetable, so that it only has one record for the corresponding cirque outline
        arcpy.Delete_management(volumetable)

        ##Calcualte the Hypsometric max and Hypsometric intergal
        array = arcpy.RasterToNumPyArray(cirqueDTM,"","","",0)
        EleArr = array[array > 0].astype(int) ##Get the elevations greater than zero
        Z_min = np.min(EleArr)
        Z_max = np.max(EleArr)
        Z_mean = np.mean(EleArr)

        Hi = (Z_mean - Z_min) / (Z_max - Z_min)
        row[12] = Hi
        
        vals,counts = np.unique(EleArr, return_counts=True)
        index = np.argmax(counts)
        hypo_max = vals[index]
        row[11] = hypo_max

        #update cursor
        cursor.updateRow(row)
        i += 1
        arcpy.AddMessage("Finished " + str(i) + " cirques")

#delete cursor variables        
del row, cursor













