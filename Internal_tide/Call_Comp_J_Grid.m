% Call the dJ_Nycander_grid function for a specified part of the grid and
% save the resulting dh and dJ 
%
clearvars; 
clc; 
close all;
% 
% Select the bounding box
lon_min = 80; %1234;
lat_min = 0;  %4321;
bbox = [lon_min lon_min + 10; lat_min lat_min + 10];

% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency
%
% Set mindepths
MinDepth = 100; % m 

% n data filename (only .mat file atm)
N_file = 'E:\Global_Data\WOD_CTD_Casts\Gridded_N_values.mat';

% the gridded bathymetry data file
bathyfile = ['E:\Global_Data\SRTM30_PLUS_w_Abyssal_Hills\' ...
             'bathy_SSG_1_120_GLOBAL_landmask.nc'];

%'bathy_SSG_1_120_GLOBAL_landmask.nc';

% save info filename
S_filename = ['IT_info_gridded_' num2str(lon_min) '_' num2str(lat_min) '.mat'];

%%%%%%%%%%Calculation region (should not need to alter below)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lon = ncread(bathyfile,'lon');
lat = ncread(bathyfile,'lat');
bathy = ncread(bathyfile,'bathy');
% Get a sufficiently large portion of the grid
I = find(lon > bbox(1,1) - 15 & lon < bbox(1,2) + 15);
J = find(lat > bbox(2,1) - 15 & lat < bbox(2,2) + 15);
lon = lon(I); lat = lat(J); bathy = bathy(I,J)';
    
%Bf = load(N_file); 
[lon, lat] = meshgrid(lon,lat);
%Compute_Nb_Nm_Gridded(lon_M,lat_M,B,zcontour,N,lon_N,lat_N)
%[~,Nm] = Compute_Nb_Nm_Gridded(lon,lat,-bathy,Bf.z,Bf.N,Bf.lon,Bf.lat);
% 
%save(S_filename,'Nm');
load(S_filename)
[dJx, dJy, lon_c, lat_c] = Calc_dJ_Nyc_Struc(lon,lat,bathy,Nm,...
                                             omega,-MinDepth,1e3,bbox,0);
                                   
save(S_filename,'dJx', 'dJy', 'dhx', 'dhy', 'lon_c', 'lat_c');     