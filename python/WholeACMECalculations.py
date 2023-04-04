#-------------------------------------------------------------------------------
# Name:        Whole ACME calculations
# Purpose:     Derive all cirque metrics based on Cirque outlines and a DEM
#
# Author: Yingkui Li
# This program derive cirque related metrics based on cirque outlines and a DEM
# The first step is to determine the cirque threshold points
# The second step is to derive length and width info, as well as the area, and parameters
# The third step is to detive the 3D statistics and hypsometric parameters
# Some of the codes are revised based on the ACME codes by Ramon Pellitero and Matteo Spagnolo 2016
# 
# Created:     08/24/2022
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
# This function derives cirque threshold points based on cirque outlines and DEM
# This function can also process the cirque outlines with and without major overlaps
#------------------------------------------------------------------------------------------------------------
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

        arcpy.CopyFeatures_management(joinedpoints, cirquePoints)
        ##delete the temp dataset
        arcpy.Delete_management (joinedpoints)

        return cirquePoints

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

        ##delete the temp dataset
        arcpy.Delete_management (sel_poly)
        arcpy.Delete_management (Singlepoint)

        return cirquePoints

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
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)
No_Overlap = arcpy.GetParameter(2)
OutputThresholds = arcpy.GetParameterAsText(3)
OutputLength = arcpy.GetParameterAsText(4)
OutputWidth = arcpy.GetParameterAsText(5)

#environments

spatialref=arcpy.Describe(InputCirques).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

arcpy.Delete_management("in_memory") ### Empty the in_memory


#a "list" where the name of fields from the attributed table are copied in
Fieldlist=[]
ListFields=arcpy.ListFields(InputCirques)
for x in ListFields:
    Fieldlist.append(x.baseName)

##Add new fields
if "L" in Fieldlist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "L", "LONG", 20)

if "W" in Fieldlist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "W", "LONG", 20)

if "L_W" in Fieldlist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "L_W", "FLOAT", 10,3)

if "Perimeter" in Fieldlist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "Perimeter", "LONG", 6)
if "Area_2D" in Fieldlist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "Area_2D", "LONG", 12)
if "Circular" in Fieldlist:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "Circular", "DOUBLE", 4, 3)

new_fields = ("Z_min","Z_max","Z_range","Z_mean","area_3D","Slope_mean", "Aspectmean", "Plan_clos", "TIN_A2D", "A3d_A2D", "hypsomax", "HI")
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(InputCirques, field, "FLOAT",6 ,2)


##Step 1: derive the threshold points
arcpy.AddMessage("Step 1: deriving threshold points based on cirque outlines and DEM")
result_cirque_thresholds = CirqueThresholds (InputCirques, InputDEM, No_Overlap)
arcpy.CopyFeatures_management(result_cirque_thresholds, OutputThresholds)
arcpy.Delete_management("in_memory") ### Empty the in_memory

##Step 2: derive 
arcpy.AddMessage("Step 2: deriving Perimeter, Area_2D, and Circularity")

fields = ("Perimeter", "Area_2D", "Circular", "SHAPE@LENGTH", "SHAPE@AREA")
with arcpy.da.UpdateCursor(InputCirques, fields) as cursor:
    for row in cursor:
       row[0]=row[3]
       row[1]=row[4]
       row[2]=row[3]/(2*math.pi*((math.sqrt(row[4]/math.pi))))
       cursor.updateRow(row)
del row, cursor

##Step 3: derive Length and width-related metrics 
arcpy.AddMessage("Step 3: deriving Length and width-related metrics")

cirque_length = arcpy.CreateFeatureclass_management("in_memory", "cirque_length", "POLYLINE")
#add ID field to attribute table of new shapefile to be matched with the FID of the original cirque polygon shapefile
arcpy.AddField_management(cirque_length, "ID_cirque", "INTEGER", 6)
cirque_width = arcpy.CreateFeatureclass_management("in_memory", "cirque_width", "POLYLINE")
#add ID field to attribute table of new shapefile to be matched with the FID of the original cirque polygon shapefile
arcpy.AddField_management(cirque_width, "ID_cirque", "INTEGER", 6)

#routine for all cirques and threshold thresh_points
with arcpy.da.UpdateCursor(InputCirques, ["L", "W", "L_W", "SHAPE@", "OID@"]) as cursor:
    for row in cursor:
        feature_class = arcpy.CopyFeatures_management(row[3], "in_memory/feature")
        thresh_p=arcpy.Intersect_analysis([OutputThresholds, feature_class], "in_memory/thres") #get the threshold inside the cirque
        n_intersections=arcpy.GetCount_management(thresh_p)#checking that there is only one intersection
        count=int(n_intersections.getOutput(0)) #checking that there is only one intersection
        if count==1:
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
        else:
            arcpy.AddMessage("cirque %r does not have a threshold or has two or more thresholds (e.g. is part of two or more nested cirques)" %(row[4])) #this relates to the cases when a cirque might not have a threshold point inside

del row, cursor

arcpy.DeleteFeatures_management(feature_class)
arcpy.DeleteFeatures_management(thresh_p)
arcpy.CopyFeatures_management(cirque_length, OutputLength)
arcpy.CopyFeatures_management(cirque_width, OutputWidth)
arcpy.Delete_management("in_memory") ### Empty the in_memory


##Step 3: derive 3d statistics and hypsometry 
arcpy.AddMessage("Step 4: derive 3d statistics and hypsometry")

fields = ("SHAPE@", "Z_min","Z_max","Z_range","Z_mean","area_3D","Slope_mean", "Aspectmean", "Plan_clos", "TIN_A2D", "A3d_A2D", "hypsomax", "HI")
volumetable = tempRs = arcpy.env.scratchFolder + "\\volumetable.txt"

with arcpy.da.UpdateCursor(InputCirques, fields) as cursor:
    i = 0
    for row in cursor:
        cirqueDTM = ExtractByMask(InputDEM, row[0]) ##shape@
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

arcpy.Delete_management("in_memory") ### Empty the in_memory












