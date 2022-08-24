#-------------------------------------------------------------------------------
# Name:        Cirque perimeter, area and circularity calculation
# Purpose:
#
# Author:      Ramon Pellitero and Matteo Spagnolo
#
# Created:     30/04/2016
# Copyright:   (c) Ramon Pellitero and Matteo Spagnolo 2016
#-------------------------------------------------------------------------------

from __future__ import division
import arcpy
from arcpy import env
import math

#inputs
cirque = arcpy.GetParameterAsText(0)

#environments
spatialref=arcpy.Describe(cirque).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

#a "list" where the name of fields from the attributed table are copied in
list=[]
ListFields=arcpy.ListFields(cirque)
for x in ListFields:
    list.append(x.baseName)

#if fields with the same name are already in the attribute table, no need to create them again (which will crash the tool)

if "Perimeter" in list:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirque, "Perimeter", "LONG", 6)

if "Area_2D" in list:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirque, "Area_2D", "LONG", 12)

if "Circular" in list:
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(cirque, "Circular", "DOUBLE", 4, 3)

#definition of fields, including the new ones
fields = ("Perimeter", "Area_2D", "Circular", "SHAPE@LENGTH", "SHAPE@AREA")

#derive metrics
with arcpy.da.UpdateCursor(cirque, fields) as cursor:
    # For each row, evaluate the Perimeter value (index position
    #  of 0), by copying the shape@length value (index position of 5)
    #
    for row in cursor:
       row[0]=row[3]
       row[1]=row[4]
       row[2]=row[3]/(2*math.pi*((math.sqrt(row[4]/math.pi))))
       cursor.updateRow(row)









