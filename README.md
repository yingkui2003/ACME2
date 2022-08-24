# Revised ACME toolbox
This revised ArcGIS toolbox derives cirque related metrics based on cirque outlines and a DEM. The major changes and added new tools include:
(1) A new tool is added to derive the cirque threshold points based on cirque outlines and the DEM.
(2) The revised tool adds a new cirque metric, the A3D/A2D ratio, because it may not be derived correctly using the Area3D and Area2D derived from the orignal ACME tool. The reason is that
the Area3D is derived from the surface volume analysis, which is based on the TIN. In some cases, the 3D area of the TIN is smaller than the Area2D derived from the polygon or raster. 
In the revised tool, the A3D/A2D ratio is directly derived using the 3D and 2D areas from the ArcGIS surface volume analysis.
(3) The revised toolbox combines the two tools of the orginal ACME for 3D statistics and Hypsometry calculation. The original ACME toolbox derived the hypsometric max and hypsometric intergal based on the contour
line interval apporach, which is very time consuming. The revised toolbox derives the hypsometric metrics based on the mode of elevation values and the Elevation-relief ratio. The revised codes is much faster
in determining the hypsometric max and hypsometric intergal values.
(4) The revised toolbox also improve the logic of the codes and fixed the memory allocation issue when processing a large amount of cirques.
(5) The revised toolbox can work in both ArcGIS and ArcGIS Pro.
