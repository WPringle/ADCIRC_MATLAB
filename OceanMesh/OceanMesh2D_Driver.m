% OceanMesh2D calling script
%        EXECUTION 
% Driver-->Prep-->distmesh2d
clearvars; close all; clc
%wrkdir = '/Users/Keith/Desktop/';
addpath('src/','-begin');
addpath('src/ann_wrapper','-begin');

%% Set parameters%%
%%-Data-%%%%%%%%%%%%%%%%%%%%%%
contourfile  = {'E:\Indian_Ocean\Shapefiles\GSHHS_shp\f\GSHHS_f_L1.shp'}; % cell of shapefile datasets in order of priority (most important first)
meshfile     = [];           %                    % Only used if contourfile is empty. This grid must contain your land boundary as nodestring.
bathyfile    = {'E:\Global_Data\SRTM90_DEM\topo15.nc'}; % cell of bathymetry datasets in order of priority (most important first)
bbox         = [131  136;  %                      % [ MinLon MaxLon ; 
                33.5  35];  %                      %   MinLat MaxLat ]             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   % This file will be used to fill in NaNs.
%%-Edge function parameters-%%                                               
min_el       = 1/480;       %                    % Minimum edgelength in degrees.
max_el       = 0.25;         %                    % Maximum edgelength in degrees. 
dist_param   = 0.20;         %                    % Distance paramater (percent that the edgelength ...
                             %                    % should change with distance, set zero to ignore).
wl_param     = 30;           %                    % Parameter in wavelength function (set zero to ignore, see ref for equation in prep...)
slope_param  = 15;           %                    % Parameter in slope function (set zero to ignore, see ref for equation in prep...)
dt           = 0;            %                    % CFL limiter: = 0 automatic, = x CFL limiter with dt = x, NB: CFL restricted to 0.5 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    % resolved by three elements degrees (set zero to ignore in prep...).
%%--For floodplain-%%%%%%%%%%%%%%
MaxDist      = 0.40;            %                 % (only used if meshfile isempty). Maximum distance in degrees that you want to trim in from the mainland boundary.
MaxEle       = 10.0;            %                 % (only used if meshfile isempty). Maximum elevation above the sea level to grid to in meters.
bounds       = [MaxEle,MaxDist];%                 % NOTE: Both maxele and maxdist must be satisfied to trim (see trim_mesh in prep).
minL         = 7.0;             %                 % Minimum length (in number of vector points) to ignore (i.e. small lakes, ponds, etc).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%-Triangle smoothing opts-%%
itmax        = 500;         %                     % Max number of iterations allowed in distmesh (generally don't change).
plot_on      = 1;           %                     % Plot (0=no plots,1=plots tri at nscreen intervals,2=additionally plots and saves plots of edge functions)
nscreen      = 1;           %                     % Output interval in number of iterations of temporary grid files and mesh quality information.
ini_p        = [];          %                     % Initial point distribution (optional, can accelerate convergence).
fix_p        = [];          %                     % Optional fix points
num_p        = 4;           %                     % Number of threads/processors to use (<=1 for serial) (generally don't use more processors than you have).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edgefx       = '';  %EdFxnIntrp                   % The edge function interp if you already did the exact same problem.
                                           
outfiname   = ['Grid_',num2str(bbox(1,1)),'_',num2str(bbox(1,2)),'_'...
                      ,num2str(bbox(2,1)),'_',num2str(bbox(2,2))]; 
% Output grid file name a "Temp_grid" is also saved during calc at nscreen
% intervals
%% Build edge function
fprintf(1,' ------------------------------------------------------->\n') ;
disp('Entering prep...')
[p,t] = OceanMesh2D_prep(meshfile,contourfile,bathyfile,bbox,min_el,max_el,...
                         dist_param,wl_param,slope_param,dt,...
                         bounds,minL,itmax,plot_on,ini_p,fix_p,nscreen,...
                         num_p,edgefx);
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
%Pool = gcp('nocreate');
%delete(Pool);

