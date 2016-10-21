% Read lithology NETCDF, map to manning's n, 
% interpolate onto mesh and output in fort.13. format
%
% Requires: readfort14.m
%
% By William Pringle Oct 19 2016

clearvars; clc; close all

% Material type reference:
% 1     Gravel and coarser                      
% 2     Sand
% 3     Silt
% 4     Clay
% 5     Calcareous ooze
% 6     Radiolarian ooze
% 7     Diatom ooze
% 8     Sponge spicules
% 9     Mixed calcareous/siliceous ooze
% 10	Shells and coral fragments
% 11	Ash and volcanic sand/gravel
% 12	Siliceous mud
% 13	Fine-grained calcareous sediment

%% Parameters to set
% Default manning's 
N_default = 0.020; 

% Mannings n table
ntable = [  0.060 ; %1  Bunya et al. 2010 (gravel pits)
            0.030 ; %2  Bunya et al. 2010 (sand bar/beach) 
            0.022 ; %3  Kerr et al. 2013 (silty)
            0.020 ; %4  Guess
            0.020 ; %5  Guess
            0.020 ; %6  Guess
            0.020 ; %7  Guess
            0.040 ; %8  Guess
            0.020 ; %9  Guess
            0.100 ; %10 Guess
            0.035 ; %11 Guess
            0.015 ; %12 Kerr et al. 2013 (muddy)
            0.025 ; %13 Guess
                  ];
              
% CPP conversion coefficents (centre of projection and projection radius)
lon0 = 75.214667 * pi/180; lat0 = -31.172085 * pi/180;
R = 6378206.4; % earth radius
              
% seabed filename
filename  = 'seabed_lithology_v1.nc';
% fort.14 name
fort14    = '../IDIOMS_v5.18.grd';
% output fort.13 name
fort13    = 'fort.13.mannings';

%% Reading NETCDF
%ncdisp(filename)
lon       = ncread(filename,'lon');
lat       = ncread(filename,'lat');
lithology = ncread(filename,'z');

%% Making manning's n and interpolant
man = zeros(size(lithology));
for n = 1:length(ntable)
    man(lithology == n) = ntable(n);
end
%man = man';
%contourf(lon,lat,man)

% Doing the CPP conversion
lon = lon * pi/180; lat = lat * pi/180; 
x = R * (lon - lon0) * cos(lat0);
y = R * lat;
% Getting the ndgrid
[lx,ly] = ndgrid(x,y);
% Gridded Interpolant
F = griddedInterpolant(lx,ly,man,'linear');

%% Read fort.14 and interpolate
[~,VX,~,~,~,~] = readfort14( fort14 ) ;

% Doing the CPP conversion
xx = VX(:,1) * pi/180; yy = VX(:,2) * pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

% Doing the interpolation
man_mesh = F(xx,yy);

%% Writing out the fort.13
% open file
fid = fopen(fort13 ,'w');
% Print header
fprintf(fid,'%s \n','Spatial attributes description') ;
fprintf(fid,'%d \n',length(VX)) ;   
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%s \n','mannings_n_at_sea_floor') ;
fprintf(fid,'%s \n','time/(length**1/3)') ;
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%f \n',N_default) ;  % some default value
fprintf(fid,'%s \n','mannings_n_at_sea_floor') ;
%number of nodes not default
MinM = min(ntable);
nnodes = length(man_mesh( man_mesh ~= N_default & man_mesh >= MinM ));
fprintf(fid,'%d \n',nnodes) ; 
% Print out list of nodes for each
for k = 1:length(VX)
    if man_mesh(k) ~= N_default && man_mesh(k) >= MinM
        fprintf(fid,'%d \t %f \n',k,man_mesh(k));
    end
end
fclose(fid);