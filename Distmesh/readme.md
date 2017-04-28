#Distmesh
- Use the files in the current directory plus the following sources to run MATLAB code:

  Distmesh source: http://persson.berkeley.edu/distmesh/

  inpoly C MEX source: www.mathworks.com/matlabcentral/fileexchange/15410-inpoly-mex-file (note that it that it does not compile with MATLAB for post 2014 versions)

  m_map: https://www.eoas.ubc.ca/~rich/map.html

  smooth2a: https://www.mathworks.com/matlabcentral/fileexchange/23287-smooth2a

- You need to have some bathymetric dataset, e.g. global SRTM30_Plus.
Ask William Pringle for details on NETCDF version. 

- You need to have made a .map file created in SMS which describes the coastlines, ocean boundaries, and islands, or you may use a .shp (shapefile) such as the GSHHS global shapefiles and then specify the bounding box. Note that this code supports dividing maps up into multiple individual ones and looping through the sub-maps inside the main .map file to make Distmesh computationally achievable for large, high resolution meshes. Furthermore .map files and .shp files can be used together.
