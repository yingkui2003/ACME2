#-------------------------------------------------------------------------------
# Name:        CirqueAxisMetrics.py
# Purpose:     Derive the axis-related cirque metrics
#
# Author: Yingkui Li
# This program derive cirque axis-related metrics based on cirque outlines, a DEM, and Length and Width features
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

#arcpy.Delete_management("in_memory") ### Empty the in_memory
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

##define normalized exp fit norm_y = a * exp(b*normx)

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

 


##Main program
# Script arguments
InputCirques  = arcpy.GetParameterAsText(0) ##Input turning points or cross sections around the outlet points
InputDEM = arcpy.GetParameterAsText(1)
InputLength = arcpy.GetParameterAsText(2)
InputWidth = arcpy.GetParameterAsText(3)

#environments

spatialref=arcpy.Describe(InputCirques).spatialReference #get spat ref from input
arcpy.env.outputCoordinateSystem = spatialref #output coordinate system is taken from spat ref
arcpy.env.overwriteOutput = True #every new created file with the same name as an already existing file will overwrite the previous file
arcpy.env.XYTolerance= "1 Meters"
arcpy.env.scratchWorkspace=arcpy.env.scratchGDB #define a default folder/database where intermediate product will be stored

arcpy.Delete_management(temp_workspace) ### Empty the in_memory

#cellsize = arcpy.GetRasterProperties_management(InputDEM,"CELLSIZEX")
#cellsize_int = int(float(cellsize.getOutput(0)))

#cirques_copy = arcpy.env.scratchGDB + "\\cirques_copy"
#arcpy.CopyFeatures_management(InputCirques, cirques_copy)

#a "list" where the name of fields from the attributed table are copied in
Fieldlist=[]
ListFields=arcpy.ListFields(InputCirques)
for x in ListFields:
    Fieldlist.append(x.baseName)

##Axis variables
new_fields = ("Axprofclos", "Axhli", "Axasp","Axgrad") ## All float variables 4
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(InputCirques, field, "DOUBLE",10, 2)

if "Axamp" in Fieldlist:  ## count = 1
    pass
else:
    #add fieds to attribute tables to be populated with calculated values
    arcpy.AddField_management(InputCirques, "Axamp", "LONG", 10)


##axis curve-fit variables 
new_fields = ("L_Exp_A","L_Exp_B","L_Exp_R2","L_Kcurv_C","L_Kcurv_R2","W_Quad_C", "W_Quad_R2") ##float variables with high digits
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(InputCirques, field, "DOUBLE",10, 4)

##axis curve-fit variables 
new_fields = ("L_NormExp_A","L_NormExp_B","L_NormExp_R2") ##float variables with high digits
for field in new_fields:
    if field in Fieldlist:
        pass
    else:
        arcpy.AddField_management(InputCirques, field, "DOUBLE",10, 4)

##Derive New cirque metrics info 
arcpy.AddMessage("Derive axis-related cirque metrics: Axprofileclos, Ax intergal, Ax Aspect, Amplitude, Axgrad, Length curve fitting parameters, and Width curve fitting parameters")

selLength = temp_workspace + "\\selLength"
selWidth = temp_workspace + "\\selWidth"

##Check the direction and flip the length from low to high elevations
arcpy.AddMessage("----derive length axis-related cirque metrics and curve fitting parameters")
Check_If_Flip_Line_Direction(InputLength, InputDEM)
arcpy.InterpolateShape_3d(InputDEM, InputLength, temp_workspace + "\\length3D", 90) ##smoothed using 3 cell-size 
arcpy.InterpolateShape_3d(InputDEM, InputWidth, temp_workspace + "\\width3D", 90) ##smoothed using 3 cell-size 

FID_list = []
HLHI_list = []
HLAsp_list = []
P_clos_list = []
Amplitude_list = []
Axgrad_list = []

L_Exp_A_list = []
L_Exp_B_list = []
L_Exp_R2_list = []

L_NormExp_A_list = []
L_NormExp_B_list = []
L_NormExp_R2_list = []

L_Kcurv_C_list = []
L_Kcurv_R2_list = []
W_Quad_C_list = []
W_Quad_R2_list = []

with arcpy.da.SearchCursor(temp_workspace + "\\length3D", ["ID_Cirque","SHAPE@", "SHAPE@Length"]) as cursor:
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
        #arcpy.AddMessage(PointZ)
        pointH = [y - min_Z for y in PointZ]
        #arcpy.AddMessage(pointH)
        #arcpy.AddMessage(LengthfromStart)
        max_H = max(pointH)
        HArr = np.array(pointH)
        norm_HArr = HArr / max_H #* 100
        #arcpy.AddMessage(norm_HArr)
        
        LenArr = np.array(LengthfromStart)
        max_len = max(LengthfromStart)
        norm_lenArr = LenArr / max_len #* 100
        #arcpy.AddMessage(norm_lenArr)
        
        #validHArr = HArr[HArr > 0]
        #validLenArr = LenArr[HArr > 0]
        validHArr = HArr[np.logical_and(HArr > 0, LenArr > 0)]
        validLenArr = LenArr[np.logical_and(HArr > 0, LenArr > 0)]


        #valid_norm_HArr = norm_HArr[HArr > 0]
        #valid_norm_lenArr = norm_lenArr[HArr > 0]
        valid_norm_HArr = norm_HArr[np.logical_and(HArr > 0, LenArr > 0)]
        valid_norm_lenArr = norm_lenArr[np.logical_and(HArr > 0, LenArr > 0)]
        
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


        #try the normalized fit 10/09/2024
        try:
            #polyfit_results = polyfit(validLenArr, [math.log(y) for y in pointH[1:]], 1)
            
            polyfit_results = polyfit(valid_norm_lenArr, np.log(valid_norm_HArr), 1)
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

        L_NormExp_A_list.append(a)
        L_NormExp_B_list.append(b)
        L_NormExp_R2_list.append(R2)

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
                         
with arcpy.da.SearchCursor(temp_workspace + "\\width3D", ["ID_Cirque","SHAPE@"]) as cursor:
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

fields = (FcID, "Axprofclos", "Axhli", "Axasp", "Axamp", "Axgrad", "L_Exp_A","L_Exp_B","L_Exp_R2","L_Kcurv_C","L_Kcurv_R2","W_Quad_C", "W_Quad_R2", "L_NormExp_A","L_NormExp_B","L_NormExp_R2")

with arcpy.da.UpdateCursor(InputCirques, fields) as cursor:
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
             
            row[13] = L_NormExp_A_list[fid]
            row[14] = L_NormExp_B_list[fid]
            row[15] = L_NormExp_R2_list[fid]
            
            #update cursor
            cursor.updateRow(row)
        except:
            arcpy.AddMessage("There is an error in the calculation. Move to the next one")
            pass
           

#delete cursor variables        
del row, cursor

arcpy.Delete_management(temp_workspace) ### Empty the in_memory












