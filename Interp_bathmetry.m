clearvars; clc; close all;
% For CPP conversion  
lon0 = 75.214893 * pi/180; lat0 = -31.179511 * pi/180;
R = 6378206.4;

%% Load mesh
load IDIOMS_v5.16.mat

% For smoothing selection
K = knnsearch(VX,[60.620161 25.303437],'k',1000);

% Convert to CPP
xx = VX(:,1) * pi/180;  yy = VX(:,2) * pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

%% Load bathymetry if interpolating
% load E:\Global_Tidal_data\GEBCO-30-arsec-data\GEBCO_IDIOMS.mat
% [lon_s,lat_s] = meshgrid(lon,lat);
% % Convert to CPP
% lon_s = lon_s * pi/180;  lat_s = lat_s * pi/180; 
% x = R * (lon_s - lon0) * cos(lat0);
% y = R * lat_s;

%K = find(B == 0);
% Now interpolate   % Put SA as bandwidth written in SMS
% B(K) = DEM2GRD(x,y,depth,xx,yy,EToV,'K',K,'SA',6549);
% B( B < -100 ) = -100;

%% Smooth the bathymetry using GUI...
B(K) = SmoothBathSelection([],xx(K),yy(K),B(K),2d3 );
% again
%B = SmoothBathSelection( EToV,xx,yy,B,1d3 );

%% Save updated mesh
save('GEBCO_on_IDIOMS_v5.16_CA.mat','VX','B');

%% Plot mesh to see new smoothed result
tri = delaunay(xx(K),yy(K));
trimesh(tri,xx(K),yy(K),B(K))
view(2)
caxis([0 1000])