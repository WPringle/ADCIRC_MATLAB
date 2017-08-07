#OceanMesh2D
1. Use the files in the current directory plus the following package to run MATLAB code:

   m_map package: https://www.eoas.ubc.ca/~rich/map.html

2. Navigate into ann_wrapper directory and run ann_class_compile. Add the ann_wrapper directory to your path. You need to have a working mex compiler. 
  
3. Make sure your wrkdir is set correctly at the top of driver. Or add permanentely through MATLAB buttons or the command window.
  
4. You need to have some bathymetric dataset, e.g. global SRTM15_PLUS available at ftp://topex.ucsd.edu/pub/srtm15_plus/.
   At the moment up to two datasets are possible, where the filenames are written into a cell type variable. A NaN in a the first dataset will be overridden by the second dataset inside the code. 

5. You need a shapefile (.shp) such as the GSHHS global shapefiles. You can include multiple shapefiles with the filenames written into the cell type variable. 
   
6. By specifying the apppropriate bounding box that you wish to mesh the program will automatically extract all it needs from the shapefile and bathy data within the bounding box for meshing.
   You may need to experiment a few times to get the exact bounding box you will like to avoid small erroneous meshed regions appearing at the edges.

7. If you want to mesh the floodplain you require an initial coastal mesh so that the program extracts the coastline from that and meshes within the floodplain matching the nodes with the coastal mesh

8. Run the program from driver. Here you set all the input parameters and file names you require.  

9. num_p specified in the driver tells MATLAB how many processors to run with for parallel function evaluations. See MATLAB settings for set limits to processors. 
