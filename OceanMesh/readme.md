#OceanMesh2D
1. Use the files in the current directory plus the following packages to run MATLAB code:

   m_map package: https://www.eoas.ubc.ca/~rich/map.html
   TopoToolbox (v2.2):  https://github.com/csdms-contrib/topotoolbox 

2. Navigate into ann_wrapper directory and run ann_class_compile. Add the ann_wrapper directory to your path. You need to have a working      mex compiler. 
   Navigate into the topotoolbox directory and run compilemexfiles.

3. Make sure your wrkdir is set correctly at the top of driver (or add permanentely through MATLAB buttons or the command window).
  
4. You need to have some bathymetric dataset, e.g. global SRTM15_PLUS available at ftp://topex.ucsd.edu/pub/srtm15_plus/.
   The files should be in either GeoTiff or ERSI ASCII text format (except for the last bathyfile in the cell  array, which is always the    global dataset, this can be netcdf). The topotoolbox can convert these files for you. NOTE: NaN in a the first dataset will be            overridden  by the global dataset inside the code. 

5. You need a shapefile (.shp) such as the GSHHS global shapefiles or a contour extracted from the bathymetric dataset using the topotoolbox GRIDObj.contour method. A contour extracted from the bathy dataset you are using is the preferred option and will give you the best quality mesh. You can include multiple shapefiles with the filenames written into the cell type variable but they may not overlap otherwise undesriable results may occur.

6. By specifying the apppropriate bounding box in the cell array "extents" that you wish to mesh, the program will automatically extract all it needs from the shapefile and bathymetric data within the bounding box for meshing. Note: the order in the cell array 'bathyfiles' and the cell array 'extents' has to correspond with each other if you are using multiple datasets.

NOTE: You can first run the scripts interactively to build the polygons encompassing the extents of the meshing region and then skip that step so you can submit it to a queue. This is done by setting PLOT_ON=0 and making the contourfile empty. This will force OceanMesh2D to load in the polygon structure .mat file from your last interactive run with those extents. This allows you to loop through multiple datasets automatically and can be very useful.   

7. If you want to mesh the floodplain, it requires an initial coastal mesh with boundaries specified (as nodestrings). The program            extracts the land boundary from the provided mesh so that nodes match with the coastal mesh exactly.
   NOTE: The floodplain requires you close your contour segements with the land in its center of the bounding box (opposite to what you      would do for the ocean). 

8. Run the program from driver. Here you set all the input parameters and file names you require.  

9. num_p specified in the driver tells MATLAB how many processors to run with for parallel function evaluations. See MATLAB settings for set limits to processors. 
