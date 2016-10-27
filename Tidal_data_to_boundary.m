%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the boundary data from the fort.14, interpolate a tidal        %
% database onto the boundary, output in fort.15 form                     %
%                                                                        %
% Requires: readfort14.m                                                 %
%                                                                        %
% Created by William Pringle Oct 20 2016 only for TPXO8 database         %
% Updated by William pringle Oct 27 2016 for FES2014 database also       % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc; 
%% Parameters to set
% CPP conversion coefficents (centre of projection and projection radius)
lon0 = 75.214667 * pi/180; lat0 = -31.172085 * pi/180;
R = 6378206.4; % earth radius
              
% Specify limits of grid
lat_min = -75; lat_max = 14;
lon_min = 20;  lon_max = 133;

% Constituents that we want in the order that we want
const = {'m2','s2','k1','o1','n2','k2','p1','q1'};

% Specify tidal database (choose TPXO8 or FES2014) and 
% directory files are in

% TPXO8
database = 'tpxo8';
database_direc = 'E:\Global_Data\TPXO8_TIDES\';

% FES2014
%database = 'FES2014';
%database_direc = 'E:\Global_Data\AVISO_TIDES\FES2014\';

% input fort.14 name
fort14    = '../IDIOMS_v5.18.grd';
% output fort.15 name
fort15    = 'fort.15.TPXOb';
%-------------------------------------------------------------------------
% Should be no reason to change below
%-------------------------------------------------------------------------
%% Evaluate tidal database file names
if strcmp(database,'tpxo8')
    tide_grid     = [database_direc 'grid_tpxo8atlas_30_v1.nc'];
    tide_prefix   = [database_direc 'hf.'];
    tide_suffix   = '_tpxo8_atlas_30c_v1.nc';
elseif strcmp(database,'FES2014')
    tide_grid     = [database_direc 'M2_FES2014b_elev.nc'];
    tide_prefix   = database_direc;
    tide_suffix   = '_FES2014b_elev.nc';
end    

%% Get boundary info
[~,VX,~,opedat,~,~] = readfort14( fort14 ) ;

b_lon = zeros(opedat.neta,1);
b_lat = zeros(opedat.neta,1);
node_num = zeros(opedat.neta,1);
ns = 1; ne = 0;
for n = 1:opedat.nope
    ne = ne + opedat.nvdll(n);
    node_num(ns:ne) = opedat.nbdv(1:opedat.nvdll(n),n);
    b_lon(ns:ne)    = VX(node_num(ns:ne),1);
    b_lat(ns:ne)    = VX(node_num(ns:ne),2);
    ns = ne + 1;
end

% Do the CPP conversion
b_lonr = b_lon * pi/180;  b_latr = b_lat * pi/180; 
b_x = R * (b_lonr - lon0) * cos(lat0);
b_y = R * b_latr;

%% Load tide data and make vectors
%ncdisp(tide_grid);
if strcmp(database,'tpxo8')
    lon = ncread(tide_grid,'lon_z');
    lat = ncread(tide_grid,'lat_z');
elseif strcmp(database,'FES2014') 
    lon = ncread(tide_grid,'lon');
    lat = ncread(tide_grid,'lat');
    lon = double(lon); lat = double(lat);
end

lat_s = length(lat);
lon_s = length(lon);
[lx,ly] = meshgrid(lon,lat);
lon_x = reshape(lx,lat_s*lon_s,1);
lat_y = reshape(ly,lat_s*lon_s,1);

% Delete uncessecary portions
% First delete by square
I = find(lon_x < lon_min | lon_x > lon_max | ...
         lat_y < lat_min | lat_y > lat_max);
lon_x(I) = []; lat_y(I) = []; 

% Delete by range search
Krs = rangesearch([lon_x,lat_y],[b_lon, b_lat],0.25); % 0.25 deg radius
LKd = 0; Kd = [];
for i = 1:length(b_lon)
   K = Krs{i};
   J = LKd + 1:LKd + length(K);
   Kd(J) = K; 
   LKd = length(Kd);
end
Kd = unique(Kd);

% The new lon and lat vectors of data
lon_x = lon_x(Kd); lat_y = lat_y(Kd);

% Do the CPP conversion
lon_x = lon_x * pi/180;  lat_y = lat_y * pi/180; 
x = R * (lon_x - lon0) * cos(lat0);
y = R * lat_y;

%% Now interpolate into grid and write out to fort.15 type file

fid = fopen(fort15,'w');

for j = 1:length(const)

    % The current consituent filename
    if strcmp(database,'tpxo8')
        tide = [tide_prefix const{j} tide_suffix];
        % For real part
        hRe = ncread(tide,'hRe');
        % reshape to vector
        Re_now = reshape(double(hRe),lat_s*lon_s,1);
        % For imaginary part
        hIm = ncread(tide,'hIm');
        % reshape to vector    
        Im_now = reshape(double(hIm),lat_s*lon_s,1);
        % Eliminate regions outside of search area and on land
        % Linear extrapolation of ocean values will be conducted where 
        % boundary nodes fall inside a land cell of the tidal data. 
        Re_now(I) = []; Re_now = Re_now(Kd); 
        K = find(Re_now == 0); Re_now(K) = []; 
        Im_now(I) = []; Im_now = Im_now(Kd); Im_now(K) = [];   
        xx = x; yy = y; xx(K) = []; yy(K) = []; 
        % Make into complex number
        Z = Re_now - Im_now*1i;
        % Do the scattered Interpolation
        F = scatteredInterpolant(xx,yy,Z,'natural');
        BZ = F(b_x,b_y);  
        BZ = BZ/1000; % mm to metres
        %
    elseif strcmp(database,'FES2014') 
        tide = [tide_prefix upper(const{j}) tide_suffix];
        % For real part
        amp_now = ncread(tide,'amplitude');
        % reshape to vector
        amp_now = reshape(double(amp_now'),lat_s*lon_s,1);
        % For imaginary part
        phs_now = ncread(tide,'phase');
        % reshape to vector    
        phs_now = reshape(double(phs_now'),lat_s*lon_s,1);
        % Eliminate regions outside of search area and on land
        % Linear extrapolation of ocean values will be conducted where 
        % boundary nodes fall inside a land cell of the tidal data. 
        amp_now(I) = []; amp_now = amp_now(Kd); 
        K = find(isnan(amp_now)); amp_now(K) = []; 
        phs_now(I) = []; phs_now = phs_now(Kd); phs_now(K) = [];   
        xx = x; yy = y; xx(K) = []; yy(K) = []; 
        % amp to m
        amp_now = amp_now/100; % cm to metres
        % phase to radians
        phs_now = phs_now * pi /180;
        % Make into complex number
        Z = amp_now.*exp(phs_now*1i);
        % Do the scattered Interpolation
        F = scatteredInterpolant(xx,yy,Z,'natural');
        BZ = F(b_x,b_y);  
    end

    % Convert real and imaginary parts to amplitude and phase
    amp_b = abs(BZ);  
    phs_b = angle(BZ)*180/pi;
    % Convert to 0 to 360;
    % phs_b that is positive 0 - 180 stays same.
    % phs_b that is negative -180 - 0 becomes 180 - 360
    phs_b(phs_b < 0) = phs_b(phs_b < 0) + 360;
    
    % Print out interpolated results
    
    % The constituent name
    fprintf(fid,'%s \n',upper(cell2mat(const(j)))) ;
    
    % Loop over the nodes of the boundaries
    for k = 1:length(amp_b)

        fprintf(fid,'%12.6f  %12.6f \t %d \t %d \t %d \n',...
                        amp_b(k),phs_b(k),j,k,node_num(k)) ;
 
    end

end

fclose(fid);
