# ACME v2 toolbox
ACME v2 is a revision and extension of the ACME toolbox (Spagnolo, M., Pellitero, R., Barr, I.D., Ely, J.C., Pellicer, X.M., Rea, B.R., 2017. ACME, a GIS tool for automated cirque metric extraction. Geomorphology 278, 280–286.) to derive cirque metrics based on cirque outlines and a DEM. The major revisons and extensions include:

(1) Expand to 48 cirque metrics, which are grouped into the dataset (projection, DEM_RES, Method), location (Lon, Lat, Easting, Northing), size (L. W, H, CS, Perimeter, A2D, A3D), shape (L_W, L_H, W_H, A3D_A2D, Circular, Slope_mean, Slope_max, Slope_min, Plan_clos, Prof_clos), aspect (Aspectmean, Asp_east, Asp_north), altitude (Z_min, Z_max, Z_mean, Z_median, Z_mid, Hypsomax, HI), axis-related (AxProfclos, Axhli, Axasp, Axamp, L_Exp_A, A_Exp_B, L_Exp_R2, L_Kcurv_C, L_Kcurv_R2, W_quad_C, W_quad_R2), and catchment-related metrics (Maxabalt, Pctabarea). 

(2) A new tool is added to autumatically derive the cirque threshold points based on cirque outlines and the DEM. Two methods are developed to derive the threshold points automatically. The first method, mainstream exit, is to derive the threshold point as the intersection point between the cirque outline and the mainstream within the cirque. The second method, threshold midpoint, is to derive the threshold point as the middle point of the cirque threshold. The users can choose one method to derive the threshod points. The users can also provide their own threshold points. ACME v2 records the threshold point method and the DEM resolution in the attribute table for the comparision purpose.
  
(3) ACME v2 combines the two tools of the orginal ACME for 3D statistics and Hypsometry calculation. The original ACME toolbox derived the hypsometric max and hypsometric intergal (HI) based on the contour
line interval apporach, which are affected by the value of contour line interval. The revised toolbox derives the hypsometric metrics based on the highest mode of elevation values and the Elevation-relief ratio 
in determining the hypsometric max and HI values, which is much faster and not affected by the contour line interval.

(4) ACME v2 revised the plan closure calculation method to be consistent with the manual method used by Dr. I.S. Evans.

(5) ACME v2 also improved the logic of the codes and fixed the memory allocation issue when processing a large amount of cirques.

Due to the end of support of ESRI for ArcGIS 10, ACME v2 is for ArcGIS Pro 2.8 or newer only. The python codes related to the tools can be run in ArcGIS 10, but the toolbox interface is newer than ArcGIS 10, so that it can not open in ArcGIS 10. Addition work is needed to recreate the tool interface in ArcGIS 10. Please contact the developer if you need to run the tools in ArcGIS 10 for more instructions.

# Major tools within the toolbox
The revised ACME toolbox includes 5 tools: 0. Derive Threshold Points from Cirque Outlines; 1. Length & Width; 2. Perimeter, Area & Circularity; 3. 3D statistics & Hypsometry, and Whole ACME calculations. 

![image](https://user-images.githubusercontent.com/24683137/186520505-b6b0468f-86bf-444c-bb33-495e14091176.png)

"Derive Threshold Points from Cirque Outlines" is designed to determine the threshold points based on cirque outlines and a DEM. The idea is to determine the threshold point of a cirque as the intersection point of the cirque outline with the major stream (highest flow accumulation). This tool can handle the cirques with or without major overlaps each other. The following is the interface of this tool:

![image](https://user-images.githubusercontent.com/24683137/186521520-4b2b6140-4541-418a-854d-6838f2196ac2.png)

"Length & Width" and "Perimeter, Area & Circularity" tools are the same with the orignal ACME other than the code revision for ArcGIS Pro.

![image](https://user-images.githubusercontent.com/24683137/186522003-4b974e5d-de8a-48be-830a-80d4e647ee9f.png)

![image](https://user-images.githubusercontent.com/24683137/186522069-cd3e6c50-17c5-49e9-acc8-414f8307d090.png)

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
