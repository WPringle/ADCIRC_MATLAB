%% OceanMesh2D calling script
%% EXECUTION Driver-->Prep-->distmesh2d
% See header of OceanMesh_Prep for explanation of inputs. 
clearvars; close all; clc
%wrkdir = '/Users/Keith/Desktop/OceanMesh_20170805/';                      %--> for myself
wrkdir=  '/home2/krober10/DistMeshRuns/OceanMesh_20170805/';               %--> for CRC
path(path,[wrkdir,'m_map']);
path(path,[wrkdir,'src']);
path(path,[wrkdir,'/src/ann_wrapper']);
addpath(genpath([wrkdir,'topotoolbox-master']));
%% Set parameters%%
contourfiles  = {'datasets/ne_se_atl.shp',...
    'datasets/se_atl_crm_v1.shp',...
    'datasets/fl_east_gom_crm_v1.shp',...
    'datasets/central_gom_crm_v1.shp',...
    'datasets/western_gom_crm_v1.shp',...
    'datasets/srtm15_ec17.shp',...
    'datasets/NOAASLR2.shp',...
    'datasets/ne_se_fl_centgom_westgom_combined.shp'};                     % Cell of shapefile datasets in consistent order with bathyfile and extents

bathyfile    = {'datasets/ne_atl_crm_v1.asc',...
    'datasets/se_atl_crm_v1.asc',...
    'datasets/fl_east_gom_crm_v1.asc',...
    'datasets/central_gom_crm_v1.asc',...
    'datasets/western_gom_crm_v1.asc',...
    'datasets/srtm15_ec17.tif',...
    'datasets/ne_comb.tif',...                                             % Cell-array of GeoTiffs, must be negative below sea level and positive above sea level.
    'datasets/topo15_compressed.nc'};                                      % This file will be used to fill in NaNs.

edgefx       =  '';%'EdFxnIntrp';                                              % The edge function interp if you already did the exact same problem.

%--> Extents are in this format-
%  [minLon maxLon
%   minLat maxLat]
extents{1}  = [-78.000  -67.5; 39 46];
extents{2}  = [-82.9996 -65.0000; 30.34 40.3];
extents{2}  = [-86 -83; 30.34 33.34];
extents{2}  = [-81.5 -79; 30.36 33.30];
extents{3}  = [-87.0000 -79.0000; 24.0000 31.0000];
extents{4}  = [-94 -87 ; 25 31];
extents{5}  = [-98 -93 ; 25 31];
extents{6}  = [-100 -15 ; 5 65];
extents{7}  = [-74.259354+1e-4, -73.800653-1e-4; 40.577485+1e-3 40.949862-1e-4];
extents{8}  = [-98.0 -67 ;25 46];

for DriverCounter = 2 : 2                                                  %-->Loop over all the regions you would like to mesh!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Building block ',num2str(DriverCounter)])
    
    %%-Edge function parameters-%%                            
    min_el       = 60/111e3;     %                                          % Minimum edgelength in degrees (if negative then use traditional linear distance function, otherwise use feature size).
    max_el       = 25e3/111e3;   %                                          % Maximum edgelength in degrees.
    max_el_ns    = 240/111e3;    %                                          % Maximum edgelength on the coastline (only relevant if min_el is positive).
    dist_param   = 0.25;         %                                          % Distance paramater (percent that the edgelength should change between neighbors).
    wl_param     = 15;            %                                          % Parameter in wavelength function (set zero to ignore, see ref for equation in prep...)
    slope_param  = 0;           %                                          % Parameter in slope function (set zero to ignore, see ref for equation in prep...)
    channel_param= 0.50e-6;      %                                          % Parameter in channel scale edge function (set zero to ignore, see ref for eq in prep...)
    dt           = 2.0;          %                                          % CFL limiter: = 0 automatic, = x CFL limiter with dt = x, NB: CFL restricted to 0.25
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%--For floodplain--%%%%%%%%%%%
    MaxDist      = 0.40;         %                                          % Maximum distance in degrees that you want to trim in from the mainland boundary.
    MaxEle       = 10.0;         %                                          % Maximum elevation above the sea level to grid to in meters.
    bounds       = [MaxEle,...   % 
        MaxDist];                %                                          % NOTE: Both maxele and maxdist must be satisfied to trim (see trim_mesh in prep).
    feature_vecs = {'datasets/se_atl_crm_v1.shp',...
        'datasets/se_atl_crm_v1_1m.shp',...
        'datasets/se_atl_crm_v1_2m.shp'};                                   % A cell array of vector features (shapefiles) that you want to constrain in the triangulation.
                                                                            % The first feature_vec is always the coastline or 0-m contour. You can run the application with a bbox and a set of feature vectors.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%-Triangle smoothing opts--%%
    itmax        = 100;         %                                           % Max number of iterations allowed in distmesh (generally don't change).
    plot_on      = 1;           %                                           % Plot (0=no plots,1=plots tri at nscreen intervals,2=additionally plots and saves plots of edge functions)
    nscreen      = 5;           %                                           % Output interval in number of iterations of temporary grid files and mesh quality information.
    ini_p        = [];          %                                           % Initial point distribution (optional, can accelerate convergence).
    fix_p        = [];          %                                           % Optional fix points
    num_p        = 2;           %                                           % Number of threads/processors to use (<=1 for serial) (generally don't use more processors than you have).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bbox         = extents{DriverCounter};
    %contourfile={contourfiles{DriverCounter}};
    contourfile='';
    
    outfiname    = ['Grid_',num2str(bbox(1,1)),'_',num2str(bbox(1,2)),'_'...% Output grid file called "Temp_grid" is also saved during calc at nscreen intervals
        ,num2str(bbox(2,1)),'_',num2str(bbox(2,2))];
    
    %% Build edge function
    fprintf(1,' ------------------------------------------------------->\n') ;
    disp('Entering prep...')
    [p,t] = OceanMesh2D_prep(contourfile,bathyfile,feature_vecs,bbox,min_el,max_el,...
        max_el_ns,dist_param,wl_param,slope_param,channel_param,dt,...
        bounds,itmax,plot_on,ini_p,fix_p,nscreen,...
        num_p,edgefx,DriverCounter,outfiname);
    disp('Exiting prep...')
    fprintf(1,' ------------------------------------------------------->\n') ;
    %% Save outputs to .mat file
    if ~isempty(p)
        fprintf(1,' ------------------------------------------------------->\n') ;
        disp(['Saving mat file of final grid with filename ',outfiname,'...']);
        save([outfiname '.mat'],'p','t');
        disp('Saving fort.14 file');
        writefort14( [outfiname '.grd'], t, p, p(:,1)*0, [] , [] ,'grid')
        disp(['Completed ',num2str(DriverCounter),' blocks.']);
    else
        fprintf(1,' ------------------------------------------------------->\n') ;
        disp('FATAL: Read output on prompt for guidance.');
    end
end

