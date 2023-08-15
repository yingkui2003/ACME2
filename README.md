# ACME v2 toolbox
ACME v2 is a revision and extension of the ACME toolbox (Spagnolo, M., Pellitero, R., Barr, I.D., Ely, J.C., Pellicer, X.M., Rea, B.R., 2017. ACME, a GIS tool for automated cirque metric extraction. Geomorphology 278, 280–286.) to derive cirque metrics based on cirque outlines and a DEM. The major revisons and extensions include:

(1) Expand to 48 cirque metrics, which are grouped into location (Lon, Lat, Easting, Northing), size (L, W, H, CS, Perimeter, A2D, A3D), shape (L_W, L_H, W_H, A3D_A2D, Circular, Slope_mean, Slope_max, Slope_min, Slpgt33, Slplt20, Slp20to33, Plan_clos, Prof_clos), aspect (Aspectmean, Asp_east, Asp_north), altitude (Z_min, Z_max, Z_mean, Z_median, Z_mid, Hypsomax, HI), axis-related (AxProfclos, Axhli, Axasp, Axamp, L_Exp_A, A_Exp_B, L_Exp_R2, L_Kcurv_C, L_Kcurv_R2, W_quad_C, W_quad_R2), and catchment-related metrics (Maxabalt, Pctabarea). ACME v2 also records 3 dataset-related attributes (projection, DEM_RES, Method), so that the derived metrics can be used for comparsion of different datasets and methods to derive the cirque threshold focus point.  

(2) A new tool is added to autumatically derive the cirque threshold focus points based on cirque outlines and the DEM. Two methods are developed to derive the threshold focus points automatically. The first method, mainstream exit, is to derive the threshold focus point as the intersection point between the cirque outline and the mainstream within the cirque. The second method, threshold midpoint, is to derive the threshold focus point as the middle point of the cirque threshold. The users can choose one method to derive the threshod points. The users can also provide their own threshold focus points. ACME v2 records the threshold focus point method and the DEM resolution in the attribute table for the comparision purpose.
  
(3) ACME v2 combines the two tools of the orginal ACME for 3D statistics and Hypsometry calculation. The original ACME toolbox derived the hypsometric max and hypsometric intergal (HI) based on the elevation bin
apporach, which are affected by the value of elevation bin. The revised toolbox derives the hypsometric metrics based on the highest mode of elevation values and the Elevation-relief ratio 
in determining the hypsometric max and HI values, which is faster and not affected by the contour line interval.

(4) ACME v2 revises the plan closure calculation method to be consistent with the manual method used by Dr. I.S. Evans.

(5) ACME v2 also improves the logic of the codes and fixed the memory allocation issue when processing a large amount of cirques.

Due to the end of support of ESRI for ArcGIS 10, ACME v2 is for ArcGIS Pro 2.8 or newer only. The python codes related to the tools can be run in ArcGIS 10, but the toolbox interface is newer than ArcGIS 10, so that it can not open in ArcGIS 10. Addition work is needed to recreate the tool interface in ArcGIS 10. Please contact the developer if you need to run the tools in ArcGIS 10 for more instructions.

# Major tools within the toolbox
ACME v2 includes two toolsets. The first one is the Step by Step tools, which include 4 tools: 0. Derive Thresholds from Cirque Outlines, 1. Record Dataset Info and Derive Location, Size, Shape, Aspect, and Altitude Metrics, 2. Derive Axis-related Metrics, and 3. Derive Catchment Metrics. The second toolset is for the Whole Process, including two tools: Whole Calculations with Auto-derived Thresholds and Whole Calculations With Specified Thresholds.

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/0187cfea-ef77-43a3-b2ed-248abd4ba4db)


"Derive Threshold Points from Cirque Outlines" derives cirque threshold focus points based on DEM and cirque polygons (outlines). The inputs include DEM, cirque outlines (polygons), and the method to derive the threshold focus points. This tool provides two methods, mainstream exit and threshold midpoint, to derive cirque threshold focus points. The output is the derived cirque threshold focus points. Make sure that the DEM and cirque outlines have the same projected coordinate system (UTM or other projected systems, not the latitudes and longitudes). The following is the interface of this tool:

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/1dd0d969-6210-4e62-8772-2251ae76faaa)




"Record Dataset Info and Derive Location, Size, Shape, Aspect, and Altitude Metrics" records dataset-related attributes and derives location, size, shape, aspect, and altitude parameters related to each cirque. The dataset-related attributes include projection, DEM_RES and the threshold focus method. The location metrics include lon, lat, easting and northing. The size parameters include L (length), W (width), H (height), CS (cirque size), Perimeter, A2D (2D area), and A3D (3D area). The shape parameters include L_W (L/W ratio), L_H (L/H ratio), W_H (W/H ratio), A3D_A2D (A3D/A2D ratio), Circular (circularity), Slope_mean, Slope_max, Slope_min, Slpgt33, Slplt20, Slp20to33, Plan_clos (plan closure), and Prof_clos (profile closure). The aspect metrics include Aspectmean, Asp_east, and Asp_north. The altitude parameters include Z_min, Z_max, Z_mean, Z_median, Z_mid (middle altitude), Hypsomax, and HI. This tool includes three inputs: Input Cirque Outlines, Input DEM, and Input Cirque Thresholds, and four outputs: Output Length Features, Output Width Features, Output Mid-Alt Contours, and Output Cirques. Make sure that the DEM, cirque outlines, and cirque thresholds have the same projected coordinate system (UTM or other projected systems, not the latitudes and longitudes).

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/206c5a13-098e-4d8a-837f-c8422546670e)


"Derive Axis-related Metrics" derives the axis-related cirque metrics, including Axprofclos, Axhli, Axasp, Axgrad and Axamp, and the axis-related curve-fitting parameters, including L_Exp_A, L_Exp_B, L_Exp_R2, L_Kcurv_C, L_Kcurv_R2, W_Quad_C, and W_Quad_R2., This tool includes two inputs: Input Cirque Outlines and Input DEM. The derived metrics are directly saved to the attribute table of the input cirque outlines. 

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/f50b9b87-7b81-4400-b4d7-4a05fe20d53f)

"Derive Catchment Metrics" derives the two metrics, Maxabalt and Pctabarea, related to the catchment area of each cirque. This tool includes two inputs: Input Cirque Outlines and Input DEM. The derived metrics are directly saved to the attribute table of the input cirque outlines.

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/12d7e4ae-e2b8-4c19-ab05-856621128621)



"Whole Calculations with Auto-derived Thresholds" derives the 48 cirque metrics (dataset, location, size, shape, aspect, altitude, axis-related and catchment-related parameters) based on the input cirque outlines and a DEM. First, the cirque threshold points are automatically derived based on the section of the two methods (mainstream exit and threshold midpoint). Then, these cirque threshold points are used to derive threshold-related cirque metrics. Note that if the users have their own threshold points, use the other tool or step-by-step tools to derive cirque metrics., This tool includes three inputs: Input Cirque Outlines, Input DEM, and the threshold point method, and five outputs: Output threshold points, Output Length Features, Output Width Features, Output Mid-Alt Contours, and Output Cirques.

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/26a0a2ee-0be5-4570-a3a7-7b42444bd9a3)

" Whole Calculations With Specified Thresholds" derives the 48 cirque metrics (dataset, location, size, shape, aspect, altitude, axis-related and catchment-related parameters) based on the input cirque outlines, DEM, and cirque threshold points. Note that the user needs to provide the cirque thresholds, which can be derived by manual digitization or step 0 in the step-by-step tools., This tool includes three inputs: Input Cirque Outlines, Input DEM, and Input threshold points, and four outputs: Output Length Features, Output Width Features, Output Mid-Alt Contours, and Output Cirques.

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/47e071a1-d740-411e-9df6-b8172d8f65aa)


# How to download and use this toolbox in ArcGIS or ArcGIS Pro
The github site includes the revised ACME toolbox (tbx) file and a python folder, including all python source codes associated with these tools. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/41865bc2-f44f-4fa9-800d-2bd5f212fe35)

A zip file of the while github folder will be downloaded to the local computer. Unzip this file will create a ACME-v2-main folder with both the tbx file and the python folder and source code files. The user can use this toolbox, check the codes, and comtinue imporving this toolbox. Note that the codes for each tool are imported into the toolbox, so that the toolbox can be run just with the tbx file only. However, the users need to export the codes first to revise the codes. 

The toolboxes and tools have been tested in from ArcGIS Pro 2.8 to 3.1. Errors may occur if using other versions of ArcGIS Pro. The python codes related to the tools can be run in ArcGIS 10, but the toolbox interface is newer than ArcGIS 10, so that it can not open in ArcGIS 10. Addition work is needed to recreate the tool interface in ArcGIS 10. 

Please report any errors or questions to Yingkui Li (yli32@utk.edu).

# Cite this work
Li et al., in preparation. TBD

# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao
