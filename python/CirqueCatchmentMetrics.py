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



##Main program
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)

#environments

spatialref=arcpy.Describe(InputCirques).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

arcpy.Delete_management("in_memory") ### Empty the in_memory

#cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
#cellsize_int = int(float(cellsize.getOutput(0)))
##Obtain the cellsize of the DEM
dem = arcpy.Raster(InputDEM)
cellsize = dem.meanCellWidth
cellsize_int = int(cellsize)

#cirques_copy = arcpy.env.scratchGDB + "\\cirques_copy"
#arcpy.CopyFeatures_management(InputCirques, cirques_copy)

#a "list" where the name of fields from the attributed table are copied in
Fieldlist=[]
ListFields=arcpy.ListFields(InputCirques)
for x in ListFields:
    Fieldlist.append(x.baseName)

##Catchment variables
if "Maxabalt" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "Maxabalt", "LONG", 10)

if "Pctabarea" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "Pctabarea", "DOUBLE",10, 2)

##Derive New cirque metrics info 
arcpy.AddMessage("Derive catchment-related cirque metrics: Maxabalt and Pctabarea")

##Derive the upstream catchment for each cirque threshold
arcpy.AddMessage("----Deriving the upstream catchment for each cirque")

fillDEM =Fill(InputDEM)  ##Fill the sink first
fdir = FlowDirection(fillDEM,"NORMAL") ##Flow direction
facc = FlowAccumulation(fdir) ##Flow accmulation

#No_Overlap = False ##Also use the one by one method to prevent the potential issues with cirque overlaps
#cirque_thresholds = CirqueThresholdsFcc (cirques_copy, facc, No_Overlap)

##Derive the streams and cross sections for cirque thresholds
Singlepoint = "in_memory\\Singlepoint"
Streams = "in_memory\\Streams"
FCselected = "in_memory\\FCselected"  
Singlestream = "in_memory\\Singlestream"
SingleWs = "in_memory\\SingleWs"
#singleBND = "in_memory\\singleBND"

#Get the minimum area of the cirque dataset
cirqueArr = arcpy.da.FeatureClassToNumPyArray(InputCirques, "SHAPE@AREA")
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

#FcID = arcpy.Describe(cirque_thresholds).OIDFieldName

#crosssections = cross_sections(cirque_thresholds, FcID, TmpStream, "MAX", 50, 120)
##Convert cirque outline polygon to polylines
cirque_lines = "in_memory\\cirque_lines"
arcpy.PolygonToLine_management(InputCirques, cirque_lines, "IGNORE_NEIGHBORS")

##Derive the watershed for cirque thresholds
TotalWS = arcpy.CreateFeatureclass_management(arcpy.env.scratchGDB, "TotalWS","POLYGON","","","",spatialref)
arcpy.AddField_management(TotalWS, "ID_cirque", "INTEGER", 6)

countResult = arcpy.GetCount_management(cirque_lines)
count = int(countResult.getOutput(0))
#FcID = arcpy.Describe(cirque_lines).OIDFieldName
#FcID = arcpy.Describe(inputFCcopy).OIDFieldName
FcID = "ORIG_FID"

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

with arcpy.da.UpdateCursor(InputCirques, fields) as cursor:
    i = 0
    for row in cursor:
        ##Get the watershed that is intersect with row[0] ##Need to check this one!! may need to do the threshold first
        #arcpy.SpatialJoin_analysis(TotalWS, row[0], selWs, "JOIN_ONE_TO_ONE", "KEEP_COMMON", '#', "INTERSECT", "1 Meters", "#")
        ##Dissolve the watershed
        #arcpy.Dissolve_management(selWs, "in_memory\\dissolve_selWs")
        #try:
        query = "ID_cirque = "+ str(row[6])
        #arcpy.AddMessage(query)
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
        #except:
        #    arcpy.AddMessage("There is an error in the calculation. Move to the next one")
        #    pass

        i += 1
        arcpy.AddMessage("Finished cirque #" + str(i))

#delete cursor variables        
del row, cursor

arcpy.Delete_management(TotalWS) 

arcpy.Delete_management("in_memory") ### Empty the in_memory











