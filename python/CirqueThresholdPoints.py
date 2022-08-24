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
import scipy
from scipy.spatial import cKDTree as KDTree
from scipy import ndimage
import arcpy.cartography as CA


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
    

def CirqueThresholds (cirque_polys, dem, no_overlap = 1):
    #get the count of autopolys
    #make copies of the data to in-memory data
    cirquepolys = "in_memory\\cirquepolys"

    cirquePoints = "in_memory\\cirquePoints"

    arcpy.CopyFeatures_management(cirque_polys, cirquepolys)

    fillDEM =Fill(dem)  ##Fill the sink first
    fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
    facc = FlowAccumulation(fdir) ##Flow accmulation
    
    countResult = arcpy.GetCount_management(cirquepolys)
    count = int(countResult.getOutput(0))
    FcID = arcpy.Describe(cirquepolys).OIDFieldName


    ##if the cirques do not overlap, do not need to loop, can just run once to extract all cirque threshold points 
    ##However, if cirques have overlaps, maybe there are some problems,
    ##Need to be tested for efficency
    ##Commented by Yingkui on 7/14/2022
    cirquePoints = arcpy.CreateFeatureclass_management("in_memory", "cirquePoints","POINT","","","",cirque_polys)

    if no_overlap:
        ##if without cirque overlaps, can do the whole analysis without looping for each cirque polygon
        ##revised and tested by Yingkui on 7/15/2022 and sucessful!!! but some times one polygon have multiple points
        outZonalStatistics = ZonalStatistics(cirquepolys, FcID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
        outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part

        arcpy.RasterToPoint_conversion(outPour, cirquePoints, "VALUE")

        ##Check if one cirque produces multiple points by the same grid_code and mark the points so that the user can choose the correct one later
        ##Get the existing fields in cirquepolys
        exist_no_use_fields = [f.name for f in arcpy.ListFields(cirquepolys)] #List of current field names in outline layer

        arcpy.AddField_management(cirquepolys, 'CirqueID', 'Long', 6) ##OFID is not FID, but it is related to FID
        arcpy.CalculateField_management(cirquepolys,"CirqueID",str("!"+str(arcpy.Describe(cirquepolys).OIDFieldName)+"!"),"PYTHON_9.3")

        joinedpoints = "in_memory\\joinedpoints"
        arcpy.SpatialJoin_analysis(cirquePoints, cirquepolys, joinedpoints, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
        
        arr = arcpy.da.FeatureClassToNumPyArray(joinedpoints, 'CirqueID')
        Cirque_IDs = np.array([item[0] for item in arr])
        unique_CirqueIDs, counts = np.unique(Cirque_IDs, return_counts=True)

        multipeIDs = unique_CirqueIDs[counts > 1] ##only find the duplicated codes
        #arcpy.AddMessage(multipeIDs)
        if len(multipeIDs) > 0: ##check if there is multiple points for some cirques
            arcpy.AddMessage("Some cirques include multiple points! These points are flagged as 1 in the attrbute table")
            arcpy.AddField_management(joinedpoints, 'Flag', 'Long', 6) ##OFID is not FID, but it is related to FID
            with arcpy.da.UpdateCursor(joinedpoints,("CirqueID","Flag")) as cursor:
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
        arcpy.DeleteField_management(joinedpoints, ["Join_Count", "TARGET_FID", "Cirque_ID", "grid_code"])

        return joinedpoints

    else:
        sel_poly = "in_memory\\sel_poly"
        Singlepoint = "in_memory\\Singlepoint"
        ##Need to loop for each cirque polygon to make sure that each cirque polygon has a threshold points
        countResult = arcpy.GetCount_management(cirquepolys)
        count = int(countResult.getOutput(0))
        FcID = arcpy.Describe(cirquepolys).OIDFieldName
        for i in range (count):
            query = FcID +" = "+ str(i+1)  ##THe FID starts from 1 for the in_memory data; Which is different from shape (FID from 0) and feature class in the geodatabase (objectID from 1)
            arcpy.AddMessage(query)
            arcpy.Select_analysis(cirquepolys, sel_poly, query)

            PolyID = arcpy.Describe(sel_poly).OIDFieldName
            
            outZonalStatistics = ZonalStatistics(sel_poly, PolyID, facc, "MAXIMUM") #Find the maximum flowaccumulation point on the moriane feature
            outPour = Con(facc == outZonalStatistics,facc)  ##Determine the highest flowaccumuation part
            
            arcpy.RasterToPoint_conversion(outPour, Singlepoint, "VALUE")
            ##Make sure only one point produced, otherwise, add a new field with flag
            pntcountResult = arcpy.GetCount_management(Singlepoint)
            pntcount = int(pntcountResult.getOutput(0))
            if pntcount > 1:
                arcpy.AddMessage("Some cirques include multiple points! These points are flagged as 1 in the attrbute table")
                arcpy.AddField_management(Singlepoint, 'Flag', 'Long', 6) ##OFID is not FID, but it is related to FID
                with arcpy.da.UpdateCursor(Singlepoint,("Flag")) as cursor:
                    for row in cursor:
                        row[1] = 1
                        cursor.updateRow(row)
                #delete cursors
                del row, cursor
            ##Append the single point to 
            arcpy.Append_management(Singlepoint, cirquePoints, "NO_TEST") 

        return cirquePoints


##--main program
##Main program
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)
No_Overlap = arcpy.GetParameter(2)
OutputCirques = arcpy.GetParameterAsText(3)

result_cirque_thresholds = CirqueThresholds (InputCirques, InputDEM, No_Overlap)

arcpy.CopyFeatures_management(result_cirque_thresholds, OutputCirques)


