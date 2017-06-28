% Distmesh calling script
clearvars; close all; clc
%% Scripts that you need
% inpoly.m (from mathworks).
% m_map toolbox.
% distmesh2d package
% limgradStruct (included in distribution)
% readfort14.m and writefort14.m   (included in distribution)
% extdom_edges.m, extdom_polygon.m (included in distribution)

% %% Setting path and compling InPolygon if required
%  wrkdir = '/Users/Keith/Desktop/';
%  path(path,[wrkdir,'m_map'])
%  path(path,[wrkdir,'DistMeshRuns/poly_stuff']);
%  path(path,[wrkdir,'distmesh']);
%  path(path,[wrkdir,'ADCIRC_MATLAB-master/ADCIRC_MATLAB-master/']);
% path(path,'/afs/crc.nd.edu/user/k/krober10/Downloads/ADCIRC_MATLAB-master/ADCIRC_MATLAB-master')
%path(path,'/afs/crc.nd.edu/user/k/krober10/Downloads/distmesh')
%path(path,'/afs/crc.nd.edu/user/k/krober10/Downloads/m_map')
% path(path,'/afs/crc.nd.edu/user/k/krober10/Downloads')
path(path,'E:\Indian_Ocean\Shapefiles\GSHHS_shp\f');
path(path,'E:\Global_Data\SRTM90_DEM');
%% Set parameters
contourfile = {'GSHHS_f_L1.shp'};                % (Only used if meshfile is empty). n-m contour where n is a number representing height...
meshfile    = '';                                % (Only used if contourfile is empty). This grid must contain your land boundary.
%bathyfile   = '/afs/crc.nd.edu/user/k/krober10/DistMeshRuns/ForTidePaper/topo15_compressed.nc';           % Bathymetry dataset (make sure to change the variable name of bathymetry in General_dismesh_fp2).
bathyfile    = 'topo15.nc';                      % Bathymetry dataset (make sure to change the variable name of bathymetry in prep).
min_el       = 1/120;                            % Minimum edgelength in degrees.
max_el       = 0.25;                             % Maximum edgelength in degrees. 
dist_param   = 0.20;                             % Distance paramater (percent that the edgelength ...
                                                 % should change with distance, set zero to ignore).
wl_param     = 200;                              % Parameter in wavelength function (set zero to ignore, see ref for equation in prep...)
slope_param  = 20;                               % Parameter in slope function (set zero to ignore, see ref for equation in prep...)
itmax        = 100;                               % Maximum number of iterations allowed in distmesh
MaxDist      = 0.40;                             % (only used if meshfile isempty). Maximum distance in degrees that you want to trim in from the mainland boundary.
MaxEle       = 10.0;                             % (only used if meshfile isempty). Maximum elevation above the sea level to grid to in meters.
bounds       = [MaxEle,MaxDist];                 % NOTE: Both maxele and maxdist must be satisfied to trim. 
minL         = 40.0;                             % (only used if meshfile is empty) Minimum length (in number of vector points) to ignore within floodplain (i.e. small lakes, ponds, etc).
plot_on      = 1;                                % Plot (1=plots triangulation at nscreen intervals,2=additionally plots and saves the edge functions, 3=debug, 0=no plots)
nscreen      = 1;                                % Output interval of temporary grid files and mesh qualityt information
ini_p        = [];                               % Initial point distribution (optional, can accelerate convergence).
fix_p        = [];                               % Fixed points (used to make the floodplain conform to the coastal mesh)... 
num_p        = 4;                               % number of processors to use (<=1 for serial)
                                                % use these if you have a
%open ocean                                     % set of vertices that you do not want to move.
% bbox        = [-50     -48                    % minLon maxLon
%                 38     40.0 ];                % minLat MaxLat
% LI area
bbox        = [-100 -80
                20 30];
% Bermuda 
% bbox        = [-67 -62 
%                 31  35];
  
outfiname   = ['Grid_',num2str(bbox(1,1)),'_',num2str(bbox(1,2)),'_'...
                      ,num2str(bbox(2,1)),'_',num2str(bbox(2,2))]; 
% Output grid file name a "Temp_grid" is also saved during calc.  

%% Build edge function
fprintf(1,' ------------------------------------------------------->\n') ;
disp('Entering prep...')
[p,t] = Distmesh2D_prep(meshfile,contourfile,bathyfile,bbox,min_el,max_el,...
                        dist_param,wl_param,slope_param,bounds,minL,itmax,...
                        plot_on,ini_p,fix_p,nscreen,num_p);
disp('Exiting prep...')
fprintf(1,' ------------------------------------------------------->\n') ;
%% Save outputs to .mat file
if ~isempty(p)
    fprintf(1,' ------------------------------------------------------->\n') ;
    disp(['Saving mat file of final grid with filename',outfiname,'...']); 
    save([outfiname '.mat'],'p','t'); 
    disp('Saving fort.14 file'); 
    writefort14( [outfiname '.grd'], t, p, p(:,1)*0, [] , [] ,'grid')
else
    fprintf(1,' ------------------------------------------------------->\n') ;
    disp('FATAL: Read output on prompt for guidance.');
end
%% Close parallel pool
Pool = gcp('nocreate');
delete(Pool);

