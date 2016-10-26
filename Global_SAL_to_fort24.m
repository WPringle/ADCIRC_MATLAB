% Interpolates the global SAL term onto the mesh and outputs a fort.24
% 
% Requires: readfort14.m
% Data required: FES2004 loads. Source at:
% ftp://ftp.legos.obs-mip.fr/pub/soa/maree/tide_model/global_solution/fes2004/
% By William Pringle Oct 20 2016

clearvars; close all; clc; 
%% Parameters to set
% CPP conversion coefficents (centre of projection and projection radius)
lon0 = 75.214667 * pi/180; lat0 = -31.172085 * pi/180;
R = 6378206.4; % earth radius
              
% Constituents that we want in the order that we want
const = {'m2','s2','k1','o1','n2','k2','p1','q1'};
frequency = [0.000140518902509,0.000145444104333,0.000072921158358,...
             0.000067597744151,0.000137879699487,0.000145842317201,...
             0.000072522945975,0.000064958541129];
alpha = {'M2 SAL','S2 SAL','K1 SAL','O1 SAL',...
         'N2 SAL','K2 SAL','P1 SAL','Q1 SAL'};

% tidal database file names
tide_grid     = 'E:\Global_Data\AVISO_TIDES\SAL\load.k1.nc';
tide_prefix   = 'E:\Global_Data\AVISO_TIDES\SAL\load.';
tide_suffix   = '.nc';
% input fort.14 name
fort14    = '../IDIOMS_v5.18.grd';
% output fort.24 name
fort24    = 'fort.24.FES2004';

%% Load tide grid data 
%ncdisp(tide_grid);
lon = ncread(tide_grid,'lon');
lat = ncread(tide_grid,'lat');

Ha = ncread(tide_grid,'Ha');

lon = lon * pi/180; lat = lat * pi/180; 
x = R * (lon - lon0) * cos(lat0);
y = R * lat;

%% Get mesh info
[~,VX,~,~,~,~] = readfort14( fort14 ) ;

% Doing the CPP conversion
xx = VX(:,1) * pi/180; yy = VX(:,2) * pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

%% Now interpolate onto grid and write out to fort.24 type file

fid = fopen(fort24,'w');

for j = 1:length(const)

    % The current consituent filename
    tide = [tide_prefix const{j} tide_suffix];
    
    % Get amp and phase
    Ha = ncread(tide,'Ha');
    Hg = ncread(tide,'Hg');
    
    Hg(Hg > 180) = Hg(Hg > 180) - 360; % move to -180 - 180
    Hg = Hg*pi/180; %radians
    
    % Convert to complex number for interpolation
    z = Ha.*exp(Hg*1i);
    
    % Do the gridded Interpolation
    F = griddedInterpolant(x,y,z,'linear');
    Z = F(xx,yy);  

    % Convert back to amp and phase
    amp = abs(Z);
    phs = angle(Z)*180/pi;
    
    % Convert to 0 to 360;
    % phs_b that is positive 0 - 180 stays same.
    % phs_b that is negative -180 - 0 becomes 180 - 360
    phs(phs < 0) = phs(phs < 0) + 360;
    
    % Print out interpolated results
    
    % The constituent details
    fprintf(fid,'%s \n',cell2mat(alpha(j))) ;
    fprintf(fid,'%17.15f \n',frequency(j)) ;
    fprintf(fid,'%d \n',1) ;  
    fprintf(fid,'%s \n',upper(cell2mat(const(j)))) ;
    
    % Loop over the nodes of the mesh
    for k = 1:length(amp)
        fprintf(fid,'%d \t %12.6f  %12.6f \n',k,amp(k),phs(k));
    end

end

fclose(fid);
