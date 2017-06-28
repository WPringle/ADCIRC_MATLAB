#Distmesh
- Use the files in the current directory plus the following sources to run MATLAB code:

  inpoly.m: www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test 

  m_map package: https://www.eoas.ubc.ca/~rich/map.html

- You need to have some bathymetric dataset, e.g. global SRTM15_PLUS available at ftp://topex.ucsd.edu/pub/srtm15_plus/

- You generally need a shapefile (.shp) such as the GSHHS global shapefiles. You can include multiple shapefiles in the cell. Then by specifying the bounding box the program will automatically extract all it needs from the shapefile and bathy data within the bounding box for meshing.

- If you want to mesh the floodplain you require an initial coastal mesh so that the program extracts the coastline from that and meshes within the floodplain matching the nodes with the coastal mesh

- Run the program from driver. Here you set all the input parameters and file names you require.  
