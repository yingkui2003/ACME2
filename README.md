# ACME v2 toolbox
ACME v2 is a revision and extension of the ACME toolbox (Spagnolo, M., Pellitero, R., Barr, I.D., Ely, J.C., Pellicer, X.M., Rea, B.R., 2017. ACME, a GIS tool for automated cirque metric extraction. Geomorphology 278, 280–286.) to derive cirque metrics based on cirque outlines and a DEM. The major revisons and extensions include:

(1) Expand to 48 cirque metrics, which are grouped into the dataset (projection, DEM_RES, Method), location (Lon, Lat, Easting, Northing), size (L. W, H, CS, Perimeter, A2D, A3D), shape (L_W, L_H, W_H, A3D_A2D, Circular, Slope_mean, Slope_max, Slope_min, Plan_clos, Prof_clos), aspect (Aspectmean, Asp_east, Asp_north), altitude (Z_min, Z_max, Z_mean, Z_median, Z_mid, Hypsomax, HI), axis-related (AxProfclos, Axhli, Axasp, Axamp, L_Exp_A, A_Exp_B, L_Exp_R2, L_Kcurv_C, L_Kcurv_R2, W_quad_C, W_quad_R2), and catchment-related metrics (Maxabalt, Pctabarea). 

(2) A new tool is added to autumatically derive the cirque threshold points based on cirque outlines and the DEM. Two methods are developed to derive the threshold points automatically. The first method, mainstream exit, is to derive the threshold point as the intersection point between the cirque outline and the mainstream within the cirque. The second method, threshold midpoint, is to derive the threshold point as the middle point of the cirque threshold. The users can choose one method to derive the threshod points. The users can also provide their own threshold points. ACME v2 records the threshold point method and the DEM resolution in the attribute table for the comparision purpose.
  
(3) ACME v2 combines the two tools of the orginal ACME for 3D statistics and Hypsometry calculation. The original ACME toolbox derived the hypsometric max and hypsometric intergal (HI) based on the contour
line interval apporach, which are affected by the value of contour line interval. The revised toolbox derives the hypsometric metrics based on the highest mode of elevation values and the Elevation-relief ratio 
in determining the hypsometric max and HI values, which is much faster and not affected by the contour line interval.

(4) ACME v2 revises the plan closure calculation method to be consistent with the manual method used by Dr. I.S. Evans.

(5) ACME v2 also improves the logic of the codes and fixed the memory allocation issue when processing a large amount of cirques.

Due to the end of support of ESRI for ArcGIS 10, ACME v2 is for ArcGIS Pro 2.8 or newer only. The python codes related to the tools can be run in ArcGIS 10, but the toolbox interface is newer than ArcGIS 10, so that it can not open in ArcGIS 10. Addition work is needed to recreate the tool interface in ArcGIS 10. Please contact the developer if you need to run the tools in ArcGIS 10 for more instructions.

# Major tools within the toolbox
ACME v2 includes two toolsets. The first one is the Step by Step tools, which include 4 tools: 0. Derive Thresholds from Cirque Outlines, 1. Derive Dataset, Location, Size, Shape, Aspect, and Altitude Metrics, 2. Derive Axis-related Metrics, and 3. Derive Catchment Metrics. The second toolset is for the Whole Process, including two tools: Whole Calculations with Auto-derived Thresholds and Whole Calculations With Specified Thresholds.

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/ca22039b-ab77-40c2-a66d-f7e2b462f436)

"Derive Threshold Points from Cirque Outlines" derives cirque threshold points based on DEM and cirque polygons (outlines). The inputs include DEM, cirque outlines (polygons), and the method to derive the threshold points. This tool provides two methods, mainstream exit and threshold midpoint, to derive cirque threshold points. The output is the derived cirque threshold points. Make sure that the DEM and cirque outlines have the same projected coordinate system (UTM or other projected systems, not the latitudes and longitudes). The following is the interface of this tool:

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/9f43d33e-0402-4487-a837-18e17f87910f)



"Derive Dataset, Location, Size, Shape, Aspect, and Altitude Metrics" derives the dataset, location, size, shape, aspect, and altitude parameters related to each cirque. The dataset metrics include projection, DEM_RES, and Method for threshold points. The location metrics include lon, lat, easting and northing. The size parameters include L (length), W (width), H (height), CS (cirque size), Perimeter, A2D (2D area), and A3D (3D area). The shape parameters include L_W (L/W ratio), L_H (L/H ratio), W_H (W/H ratio), A3D_A2D (A3D/A2D ratio), Circular (circularity), Slope_mean, Slope_max, Slope_min, Plan_clos (plan closure), and Prof_clos (profile closure). The aspect metrics include Aspectmean, Asp_east, and Asp_north. The altitude parameters include Z_min, Z_max, Z_mean, Z_median, Z_mid (middle altitude), Hypsomax, and HI. This tool includes four inputs: Input Cirque Outlines, Input DEM, Input Cirque Threshold points, and Threshold point method, and four outputs: Output Length Features, Output Width Features, Output Mid-Alt Contours, and Output Cirques. Make sure that the DEM, cirque outlines, and cirque thresholds have the same projected coordinate system (UTM or other projected systems, not the latitudes and longitudes).

![image](https://github.com/yingkui2003/ACME-v2/assets/24683137/1ace55de-07bc-4e54-89d1-302fd95a4415)




"3D statistics & Hypsometry" is the revised tool to combine the orginal ACME toolbox for 3D statistics and Hypsometry calculation. 

![image](https://user-images.githubusercontent.com/24683137/186522395-ec91c9d0-b591-4101-9c71-57a39afecbe3.png)

"Whole ACME calculations" is developed to derive all cirque metrics based on the inputs of cirque outlines and the DEM. In addtion to the cirque metrics, this tool also generate the cirque threshold points, cirque length and width features.

![image](https://user-images.githubusercontent.com/24683137/186522720-b48b6f07-6f23-4b9e-9fa2-3147f724a7bb.png)


# How to download and use this toolbox in ArcGIS or ArcGIS Pro
The github site includes the revised ACME toolbox (tbx) file and a python folder, including all python source codes associated with these tools. The user can click "Code" (green color) on the right side of the github page and choose Download Zip.

![image](https://user-images.githubusercontent.com/24683137/186519537-8c7455ae-f69d-4d26-9ef1-13965c453a92.png)

A zip file of the while github folder will be downloaded to the local computer. Unzip this file will create a RevisedACME-main folder with both the tbx file and the python folder and source code files. The user can use this toolbox, check the codes, and comtinue imporving this toolbox. Note that the codes for each tool are imported into the toolbox, so that the toolbox can be run just with the tbx file only. However, the users need to export the codes first to revise the codes. 

The toolboxes and tools have been tested in ArcGIS 10.7, 10.8, 10.9 and ArcGIS Pro 2.8 and 2.9. Errors may occur if using other versions of ArcGIS or ArcGIS Pro. 

Please report any errors or questions to Yingkui Li (yli32@utk.edu).

# Cite this work
Li Y., in preparation. TBD

# Contact info
Yingkui Li

Department of Geography & Sustainability

University of Tennessee

Knoxville, TN 37996

Email: yli32@utk.edu

Website: https://geography.utk.edu/about-us/faculty/dr-yingkui-li/

Google Scholar: https://scholar.google.com/citations?user=JoNuyCMAAAAJ&hl=en&oi=ao
