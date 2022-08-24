#-------------------------------------------------------------------------------
# Name:        Cirque lenght and width calculation
# Purpose:
#
# Author:      Ramon Pellitero and Matteo Spagnolo
#
# Created:     30/04/2016
# Copyright:   (c) Ramon Pellitero and Matteo Spagnolo 2016
#
# Revised by Yingkui Li on 8/24/2022, so that the codes can work on both ArcMap and ArcGIS Pro.
#
#-------------------------------------------------------------------------------
from __future__ import division
import arcpy, arcinfo
from arcpy import env
from arcpy.sa import *
import math


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
    
#inputs
threshold_thresh_point = arcpy.GetParameterAsText(0)
cirque = arcpy.GetParameterAsText(1)
L_W_shp_folder = arcpy.GetParameterAsText(2)

#environments
spatialref=arcpy.Describe(cirque).spatialReference
arcpy.env.outputCoordinateSystem = spatialref
arcpy.env.overwriteOutput = True
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB

#a "list" where the name of fields from the attributed table are copied in
list=[]
ListFields=arcpy.ListFields(cirque)
for x in ListFields:
    list.append(x.baseName)

if "L" in list:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirque, "L", "LONG", 20)

if "W" in list:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirque, "W", "LONG", 20)

if "L_W" in list:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirque, "L_W", "FLOAT", 10,3)

cirque_length_shp = arcpy.CreateFeatureclass_management(L_W_shp_folder, "cirque_length.shp", "POLYLINE")
#add ID field to attribute table of new shapefile to be matched with the FID of the original cirque polygon shapefile
arcpy.AddField_management(cirque_length_shp, "ID_cirque", "INTEGER", 6)
cirque_width_shp = arcpy.CreateFeatureclass_management(L_W_shp_folder, "cirque_width.shp", "POLYLINE")
#add ID field to attribute table of new shapefile to be matched with the FID of the original cirque polygon shapefile
arcpy.AddField_management(cirque_width_shp, "ID_cirque", "INTEGER", 6)

#routine for all cirques and threshold thresh_points
with arcpy.da.UpdateCursor(cirque, ["L", "W", "L_W", "SHAPE@", "OID@"]) as rows:
    for row in rows:
        feature_class = arcpy.CopyFeatures_management(row[3], "in_memory/feature")
        thresh_p=arcpy.Intersect_analysis([threshold_thresh_point, feature_class], "in_memory/thres") #get the threshold inside the cirque
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

            cursor = arcpy.da.InsertCursor(cirque_length_shp, ["SHAPE@", "ID_cirque"])
            cursor.insertRow([LENGTH, row[4]])
            del cursor

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
            cursor = arcpy.da.InsertCursor(cirque_width_shp, ["SHAPE@", "ID_cirque"])
            cursor.insertRow([WIDTH, row[4]])
            del cursor

            #populate field in cirque shapefile
            row[1]=int(WIDTH.length)

            #calculate L/W
            row[2]=row[0]/row[1]

            #update rows content so that it is permanently stored
            rows.updateRow(row)


        else:
            arcpy.AddMessage("cirque %r does not have a threshold or has two or more thresholds (e.g. is part of two or more nested cirques)" %(row[4])) #this relates to the cases when a cirque might not have a threshold point inside
        arcpy.DeleteFeatures_management(feature_class)
        arcpy.DeleteFeatures_management(thresh_p)
        arcpy.Delete_management("in_memory")
















