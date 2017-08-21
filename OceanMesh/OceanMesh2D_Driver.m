% OceanMesh2D calling script
%        EXECUTION
% Driver-->Prep-->distmesh2d
clearvars; close all; clc
wrkdir = '/Users/Keith/Desktop/OceanMesh_20170805/';
path(path,[wrkdir,'m_map']);
path(path,[wrkdir,'src']);
path(path,[wrkdir,'/src/ann_wrapper']);
addpath(genpath([wrkdir,'topotoolbox-master']));
global DriverCounter
%% Set parameters%%
contourfiles  = {'datasets/ne_se_atl.shp',...
    'datasets/se_atl_crm_v1.shp',...
    'datasets/fl_east_gom_crm_v1.shp',...
    'datasets/central_gom_crm_v1.shp',...
    'datasets/western_gom_crm_v1.shp'};                         % Cell of shapefile datasets in consistent order with bathyfile and extents
meshfile     = [];                                                          % This grid must contain your land boundary as nodestring.
bathyfile    = {'datasets/ne_atl_crm_v1.asc',...
    'datasets/se_atl_crm_v1.asc',...
    'datasets/fl_east_gom_crm_v1.asc',...
    'datasets/central_gom_crm_v1.asc',...
    'datasets/western_gom_crm_v1.asc',...                        % Cell of bathymetry datasets in order of priority (most important first).
    '/Users/Keith/Desktop/crm_coasltine/topo15_compressed.nc'}; % This file will be used to fill in NaNs.
% [minLon maxLon
%  minLat maxLat]
extents{1}  = [-78.000  -76.5; 39 41];
extents{2}  = [-82.9996 -65.0000; 30.34 40.3];
extents{3}  = [-87.0000 -79.0000; 24.0000 31.0000];
extents{4}  = [-94 -87 ; 25 31];
extents{5}  = [-98 -93 ; 25 31];
%bbox          =[-89.99 -89.91;30 30.1];


for DriverCounter = 1 : length(extents)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['Building block ',num2str(DriverCounter)])
    %%-Edge function parameters-%%
    min_el       = 3/1200/4;     %                    % Minimum edgelength in degrees.
    max_el       = 0.25;         %                    % Maximum edgelength in degrees.
    dist_param   = 0.20;         %                    % Distance paramater (percent that the edgelength ...
    %                            %                    % should change with distance, set zero to ignore).
    wl_param     = 15;           %                    % Parameter in wavelength function (set zero to ignore, see ref for equation in prep...)
    slope_param  = 10;           %                    % Parameter in slope function (set zero to ignore, see ref for equation in prep...)
    confluence_scale = 1;        %                    % Use confluence scale edge function (1) or don't use it (0).
    dt           = 4;            %                    % CFL limiter: = 0 automatic, = x CFL limiter with dt = x, NB: CFL restricted to 0.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    % resolved by three elements degrees (set zero to ignore in prep...).
    %%--For floodplain--%%%%%%%%%%%%%%
    MaxDist      = 0.40;            %                 % (only used if meshfile isempty). Maximum distance in degrees that you want to trim in from the mainland boundary.
    MaxEle       = 10.0;            %                 % (only used if meshfile isempty). Maximum elevation above the sea level to grid to in meters.
    bounds       = [MaxEle,MaxDist];%                 % NOTE: Both maxele and maxdist must be satisfied to trim (see trim_mesh in prep).
    minL         = 7.0;             %                 % Minimum length (in number of vector points) to ignore (i.e. small lakes, ponds, etc).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%-Triangle smoothing opts--%%
    itmax        = 500;         %                     % Max number of iterations allowed in distmesh (generally don't change).
    plot_on      = 1;           %                     % Plot (0=no plots,1=plots tri at nscreen intervals,2=additionally plots and saves plots of edge functions)
    nscreen      = 5;           %                     % Output interval in number of iterations of temporary grid files and mesh quality information.
    ini_p        = [];          %                     % Initial point distribution (optional, can accelerate convergence).
    fix_p        = [];          %                     % Optional fix points
    num_p        = 2;           %                     % Number of threads/processors to use (<=1 for serial) (generally don't use more processors than you have).
                                                      % Some versions of MATLAB require you set this in parallel preferences.                                             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    edgefx       = '';  %EdFxnIntrp                   % The edge function interp if you already did the exact same problem.
    
    bbox         = extents{DriverCounter};
    contourfile  = {contourfiles{DriverCounter}};
    
    outfiname   = ['Grid_',num2str(bbox(1,1)),'_',num2str(bbox(1,2)),'_'...
        ,num2str(bbox(2,1)),'_',num2str(bbox(2,2))];
    % Output grid file name a "Temp_grid" is also saved during calc at nscreen
    % intervals
    %% Build edge function
    fprintf(1,' ------------------------------------------------------->\n') ;
    disp('Entering prep...')
    [p,t] = OceanMesh2D_prep(meshfile,contourfile,bathyfile,bbox,min_el,max_el,...
        dist_param,wl_param,slope_param,confluence_scale,dt,...
        bounds,minL,itmax,plot_on,ini_p,fix_p,nscreen,...
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
        close all;        
    else
        fprintf(1,' ------------------------------------------------------->\n') ;
        disp('FATAL: Read output on prompt for guidance.');
    end
end
%% Close parallel pool
Pool = gcp('nocreate');
delete(Pool);

