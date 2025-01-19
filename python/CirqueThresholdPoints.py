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

temp_workspace = "in_memory"  
if ArcGISPro:
    temp_workspace = "memory"

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
    #arcpy.AddMessage(lengths)
    if len(lengths) == 1:
        return line
    startpnts = temp_workspace + "\\startpnts"
    endpnts = temp_workspace + "\\endpnts"
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

    ##Find the shortest distance inorder until pntCount - 1
    ##sort the distlist from small to large
    distArr = np.array(distlist)
    indexArr = np.argsort(distArr)

    connectionline = arcpy.CreateFeatureclass_management(temp_workspace, "connectionline","POLYLINE","","","",line)
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
        minlength = min(start_line_length, end_line_length)
        #arcpy.AddMessage("minlength: " + str(minlength))

        mingap = min(minlength, 60)
        #arcpy.AddMessage("mingap: " + str(mingap))

        if distArr[index] < mingap:
            #arcpy.AddMessage("Connect two lines")
            array = arcpy.Array([arcpy.Point(startpntX[start],startpntY[start]),arcpy.Point(endpntX[end], endpntY[end])])
            polyline = arcpy.Polyline(array)
            new_line_cursor.insertRow([polyline])

    del new_line_cursor

    arcpy.Append_management(connectionline, line, "NO_TEST")
    arcpy.Dissolve_management(line, temp_workspace + "\\line_dissolve", "", "", "SINGLE_PART")
    ##check to delete unconnected lines
    lineArray = arcpy.da.FeatureClassToNumPyArray(temp_workspace + "\\line_dissolve",["SHAPE@LENGTH"])
    if len(lineArray) > 1:
        #arcpy.AddMessage("Delete unconnected small line sections")
        lengths = np.array([item[0] for item in lineArray])
        #arcpy.AddMessage(lengths)
        max_length = max(lengths)
        with arcpy.da.UpdateCursor(temp_workspace + "\\line_dissolve",("SHAPE@LENGTH")) as cursor:
            for row in cursor:
                if row[0] < max_length:
                    cursor.deleteRow()
        #delete cursors
        del row, cursor
    
    arcpy.CopyFeatures_management(temp_workspace + "\\line_dissolve", line)

    arcpy.Delete_management (temp_workspace + "\\line_dissolve")
    arcpy.Delete_management (connectionline)
    
    return line

def CirqueThresholds (cirque_polys, dem): ##, overlap = False):
    #get the count of autopolys
    #make copies of the data to in-memory data
    cirquepolys = temp_workspace + "\\cirquepolys"

    cirquePoints = temp_workspace + "\\cirquePoints"

    arcpy.CopyFeatures_management(cirque_polys, cirquepolys)

    ##Check if the cirques have overlaps
    '''
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
    '''
    fillDEM =Fill(dem)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
    facc = FlowAccumulation(fdir) ##Flow accmulation
    
    #countResult = arcpy.GetCount_management(cirquepolys)
    #count = int(countResult.getOutput(0))
    FcID = arcpy.Describe(cirquepolys).OIDFieldName

    b_Overlap = True ##always use the overlapped apporach
    
    cirquePoints = arcpy.CreateFeatureclass_management(temp_workspace, "cirquePoints","POINT","","","",cirque_polys)

    if b_Overlap == False:
        ##if without cirque overlaps, can do the whole analysis without looping for each cirque polygon
        ##revised and tested by Yingkui on 7/15/2022 and sucessful!!! but some times one polygon have multiple points
        outZonalStatistics = ZonalStatistics(cirquepolys, FcID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part

        arcpy.RasterToPoint_conversion(outPour, cirquePoints, "VALUE")

        ##Check if one cirque produces multiple points by the same grid_code and mark the points so that the user can choose the correct one later
        ##Get the existing fields in cirquepolys
        exist_no_use_fields = [f.name for f in arcpy.ListFields(cirquepolys)] #List of current field names in outline layer

        arcpy.AddField_management(cirquepolys, 'ID_Cirque', 'Long', 6) ##OFID is not FID, but it is related to FID
        arcpy.CalculateField_management(cirquepolys,"ID_Cirque",str("!"+str(arcpy.Describe(cirquepolys).OIDFieldName)+"!"),"PYTHON_9.3")

        joinedpoints = temp_workspace + "\\joinedpoints"
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

        arcpy.CopyFeatures_management(joinedpoints, cirquePoints)
        ##delete the temp dataset
        arcpy.Delete_management (joinedpoints)

        return cirquePoints

    else:
        sel_poly = temp_workspace + "\\sel_poly"
        sel_outline = temp_workspace + "\\sel_outline"
        outline_buf = temp_workspace + "\\outline_buf"
        Singlepoint = temp_workspace + "\\Singlepoint"
        Singlepoint_ele = temp_workspace + "\\Singlepoint_Ele"
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

            ##convert the polygon to polyline
            arcpy.PolygonToLine_management(sel_poly, sel_outline)
            #arcpy.analysis.Buffer(sel_outline, outline_buf, "10 Meters")

            PolyID = arcpy.Describe(sel_poly).OIDFieldName
            
            outZonalStatistics = ZonalStatistics(sel_poly, PolyID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature

            #LineID = arcpy.Describe(sel_outline).OIDFieldName
            #outZonalStatistics = ZonalStatistics(sel_outline, LineID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature

            outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part

            ##Snap the outPour to the outline boundary
            
            arcpy.RasterToPoint_conversion(outPour, Singlepoint, "VALUE")
            arcpy.AddField_management(Singlepoint, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID
            with arcpy.da.UpdateCursor(Singlepoint,("ID_cirque")) as cursor:
                for row in cursor:
                    row[0] = fid
                    cursor.updateRow(row)
            #delete cursors
            del row, cursor

            snapEnv1 = [sel_outline, "VERTEX", "30 Meter"]
            snapEnv2 = [sel_outline, "EDGE", "30 Meter"]
            arcpy.Snap_edit(Singlepoint, [snapEnv1, snapEnv2])           
            ##Make sure only one point produced, otherwise, add a new field with flag
            pntcountResult = arcpy.GetCount_management(Singlepoint)
            pntcount = int(pntcountResult.getOutput(0))
            if pntcount > 1:
                arcpy.AddMessage("Some cirques include multiple points! choose the lowest one as the threshold point")
                #arcpy.AddField_management(Singlepoint, 'Flag', 'Long', 6) ##OFID is not FID, but it is related to FID
                #with arcpy.da.UpdateCursor(Singlepoint,("Flag")) as cursor:
                #arcpy.CopyFeatures_management(Singlepoint, "d:\\temp\\singlepoint.shp")
                ExtractValuesToPoints(Singlepoint, dem, Singlepoint_ele, "INTERPOLATE", "VALUE_ONLY")
                pntarray = arcpy.da.FeatureClassToNumPyArray(Singlepoint_ele,"RASTERVALU")
                pntElevs = np.array([item[0] for item in pntarray])
                #arcpy.AddMessage(pntElevs)
                minEleV = int(min(pntElevs))
                with arcpy.da.UpdateCursor(Singlepoint_ele,("RASTERVALU")) as cursor:
                    i = 0
                    for row in cursor:
                        if (int(row[0]) > minEleV) or (i > 0): ##just keep one point
                            cursor.deleteRow()
                        else:
                            i += 1
                            
                #delete cursors
                del row, cursor
                ##Snap the point to the outline
                #arcpy.Snap_edit(Singlepoint_ele, [snapEnv1, snapEnv2])
                ##Append the single point to 
                arcpy.Append_management(Singlepoint_ele, cirquePoints, "NO_TEST") 
            else:
                ##Snap the point to the outline
                #arcpy.Snap_edit(Singlepoint, [snapEnv1, snapEnv2])
                arcpy.Append_management(Singlepoint, cirquePoints, "NO_TEST") 
        ##delete the temp dataset
        arcpy.Delete_management (sel_poly)
        arcpy.Delete_management (Singlepoint)

        return cirquePoints


def CirqueThresholds_midpoints(InputCirques, InputDEM):####, percentile):
    ##This method determines the cirque threshold mid-points by the lower 10th percentile of the elevation of the cirque outlines
    sel_poly = temp_workspace + "\\sel_poly"
    cirque_line = temp_workspace + "\\cirque_line"
    tmpbuf = temp_workspace + "\\tmpbuf"
    points = temp_workspace + "\\points"
    points_with_Z = temp_workspace + "\\points_with_Z"
    sel_points = temp_workspace + "\\sel_points"
    tmpline = temp_workspace + "\\tmpline"
    tmpline3D = temp_workspace + "\\tmpline3D"
    midpoint = temp_workspace + "\\midpoint"

    ##Obtain the cellsize of the DEM
    dem = arcpy.Raster(InputDEM)
    cellsize = dem.meanCellWidth
    #arcpy.AddMessage("Cell size: " + str(cellsize))

    cirquePoints = arcpy.CreateFeatureclass_management(temp_workspace, "cirquePoints","POINT","","","",InputCirques)
    arcpy.AddField_management(cirquePoints, 'ID_cirque', 'Long', 6) ##OFID is not FID, but it is related to FID

    cirqueThresholds = arcpy.CreateFeatureclass_management(temp_workspace, "cirqueThresholds","POLYLINE","","","",InputCirques)
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

        midElev = (np.max(Elevs) + np.min(Elevs)) / 2
        lowhalfElevs = Elevs[Elevs < midElev]
        #binvalue = int((maxElev - minElev) /5)
        #arcpy.AddMessage(lowhalfElevs)
        bins_num = int((midElev -np.min(Elevs))/5)

        prob, ele = np.histogram(lowhalfElevs, bins = bins_num)
        peak_prob_elev = ele[np.where(prob == prob.max())][0]

        cutoff_elev = peak_prob_elev + 2 * 5 ##add one more bin elevation (5 meter) to the peak elevation to include all potential cirque threshold

        ##Method 1: use the raster functions to get the cirque outlines within the lower elevations
        outCon = Con(cirqueDTM < cutoff_elev, 1)
        OutBndCln = BoundaryClean(outCon)
        arcpy.conversion.RasterToPolygon(OutBndCln, temp_workspace + "\\OutBndCln_poly")
        arcpy.analysis.Clip(cirque_line, temp_workspace + "\\OutBndCln_poly", tmpline)
        arcpy.MultipartToSinglepart_management(tmpline, temp_workspace + "\\tmpline_singlePart")
        arcpy.Dissolve_management(temp_workspace + "\\tmpline_singlePart", tmpline, "", "", "SINGLE_PART")

        ##test if there are more line sections
        lineArray = arcpy.da.FeatureClassToNumPyArray(tmpline,["SHAPE@LENGTH"])
        if len(lineArray) == 0: ##if no line is created
            cutoff_elev = np.percentile(Elevs, 25) 
            outCon = Con(cirqueDTM < cutoff_elev, 1)
            OutBndCln = BoundaryClean(outCon)
            arcpy.conversion.RasterToPolygon(OutBndCln, temp_workspace + "\\OutBndCln_poly")
            arcpy.analysis.Clip(cirque_line, temp_workspace + "\\OutBndCln_poly", tmpline)
            arcpy.MultipartToSinglepart_management(tmpline, temp_workspace + "\\tmpline_singlePart")
            arcpy.Dissolve_management(temp_workspace + "\\tmpline_singlePart", tmpline, "", "", "SINGLE_PART")
            
        if len(lineArray) > 1:
            #arcpy.AddMessage("multiple lines created, need to connect the lines")
            ##only connected major lines
            ##remove the small lines
            connect_major_line_sections(tmpline)

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
    arcpy.Delete_management (temp_workspace + "\\OutBndCln_poly")
    arcpy.Delete_management (temp_workspace + "\\OutBndCln_poly")
    arcpy.Delete_management (temp_workspace + "\\tmpline_singlePart")
    
    return cirquePoints, cirqueThresholds


##--main program
##Main program
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)
method = arcpy.GetParameterAsText(2)
OutputCirques = arcpy.GetParameterAsText(3)

arcpy.Delete_management(temp_workspace) ### Empty the in_memory

##make sure the projection of the input datasets is the same projected projection 
spatial_ref_cirques = arcpy.Describe(InputCirques).spatialReference
spatial_ref_dem = arcpy.Describe(InputDEM).spatialReference


#if "UTM" in spatial_ref_dem.name:
if spatial_ref_dem.linearUnitName == "Meter":
    arcpy.AddMessage("The DEM projection is: " + spatial_ref_dem.name)
else:
    arcpy.AddMessage("The unit of the DEM projection is not in meter. Please re-project the DEM to a projected coordinate system for the analysis!")
    exit()   

#if "UTM" in spatial_ref_crosssections.name:
if spatial_ref_cirques.linearUnitName == "Meter":
    arcpy.AddMessage("The cirque outline projection is: " + spatial_ref_cirques.name)
else:
    arcpy.AddMessage("The unit of the cirque outline projection is not in meter. Please re-project it to a projected coordinate system for the analysis!")
    exit()   

if spatial_ref_dem.name == spatial_ref_cirques.name:
    arcpy.AddMessage("Both DEM and cirque outlines have the same projected coordinate system: " + spatial_ref_dem.name)
else:
    arcpy.AddMessage("The DEM and cirque outlines have different map projections. Please re-project the datasets to the same projection!")
    exit()   

#arcpy.AddMessage(method)

if method == "Mainstream exit":
    result_threshold_midpoints = CirqueThresholds (InputCirques, InputDEM)
else:
    result_threshold_midpoints, cirqueThresholds = CirqueThresholds_midpoints(InputCirques, InputDEM)###, 10)
    
arcpy.CopyFeatures_management(result_threshold_midpoints, OutputCirques)

arcpy.Delete_management(temp_workspace) ### Empty the in_memory




