%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Title:       Make_IT_Fric                                      %
%  Description: Read constant contours of N values over depth     %
%               and bathymetry data to compute internal tide      %
%               friction coefficients on unstructured mesh        %
%  Inputs:      1) .mat files of N values at constant contours    %
%               2) Bathymetry data (e.g. gridded GEBCO)           %
%               3) Unstructured grid nodes                        %                           
%  Outputs:     A fort.13 formatted file for use in ADCIRC/SMS    %
%  Project:     Indian Ocean and Marginal Seas                    %
%  Author:      William Pringle                                   %
%  Created:     Oct 5 2016                                        %
%  Updated:     Oct 24 2016
%  Requires:    function - readfort14_nodes, Compute_Nb_Nm        %
%               & ADCIRC_Bath_Slope                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clearvars; 
clc; 
close all;
%
%% Set parameters and filenames here
% CPP conversion coefficents (centre of projection and projection radius)
lon0 = 75.214667 * pi/180; lat0 = -31.172085 * pi/180;
R = 6378206.4;

% Coriolis coefficient
psi = 2*7.29212d-5;

% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency

% Choose the dimensionless coefficient for the internal tide friction
C_it = 2/3;

% Minimum and maximum allowable IT_Fric 
MinFoldingTime = 6;  % hours
MaxFoldingTime = 30; % days
MinDepth       = 100; % m 

MaxIT = 1/(MinFoldingTime*3600);
MinIT = 1/(MaxFoldingTime*24*3600);

% Bathmetry data filename (only .mat file atm) 
bath_filename = 'E:\Global_Data\GEBCO-30-arsec-data\GEBCO_IDIOMS.mat';

% ADCIRC mesh filame (fort.14)
adcirc_filename = 'IDIOMS_v5.18.grd';

% n data filename (only .mat file atm)
N_filename = 'E:\Global_Data\WOD_CTD_Casts\Indian_Ocean_N_100m_processed.mat';
%%
%%%%%%%%%%Calculation region (should not need to alter below)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Load GEBCO data and interpolate grad^2 onto mesh;
% load(bath_filename); 
% 
% % CPP conversion to x, y and calculation of step size
% lonr = lon' * pi/180; latr = lat' * pi/180; 
% x = R * (lonr - lon0) * cos(lat0);
% y = R * latr;
% dx = x(2)-x(1);
% dy = y(2)-y(1);
% 
% % Evaluate slope
% [Hx, Hy] = gradient(depth,dx,dy);
% 
% % Grad^2
% H2 = Hx.^2 + Hy.^2;
% 
% % Make gridded interpolant
% [xv, yv] = ndgrid(x,y);
% F = griddedInterpolant(xv,yv,H2','linear');

%% Load ADCIRC mesh and interpolate onto
%load(adcirc_filename)
[EToV,VX,B,~,~,~] = readfort14(adcirc_filename);

xx = VX(:,1) * pi/180;  yy = VX(:,2)* pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

% Do the interpolation onto the ADCIRC mesh
[Hx,Hy] = ADCIRC_Bath_Slope( EToV,xx,yy, B ,6549 );
H2_mesh = Hx.^2 + Hy.^2;

% H2_mesh = F(xx,yy);
%
%% Load the constant contours of N values and compute Nb and Nmean
load(N_filename); 
%
[Nb,Nm] = Compute_Nb_Nm(VX(:,1),VX(:,2),B,zcontour,...
                        N,lon,lat,lon0,lat0);
                    
%% Calculate F_it from Nb, Nm, and H2_mesh
F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2))...
       .*H2_mesh/omega;
F_it = real(F_it); % in case becomes complex due to minus square root

% Compute criticality
f = psi*sind(VX(:,2));
alpha  = (omega^2 - f.^2)./(Nm.^2 - omega^2); % Nm or Nb??
gamma2 = max(H2_mesh./alpha,1);

% Normalise F_it by criticality
F_it = F_it./gamma2;

% 
%% Write out the fort.13 file
filename = 'fort.13.f_it_crit';
fid = fopen(filename,'w');
% Print header
fprintf(fid,'%s \n','Spatial attributes description') ;
fprintf(fid,'%d \n',length(B)) ;   
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%s \n','internal_tide_friction') ;
fprintf(fid,'%s \n','1/time') ;
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%f \n',0.0) ;  
fprintf(fid,'%s \n','internal_tide_friction') ;
% 
% Clipping F_it for small and large e-folding times, and small depths
%F_it(F_it > MaxIT) = MaxIT;
F_it(F_it < MinIT) = 0;
F_it(B < MinDepth) = 0;
%
% Number of nodes not default (0.0)
ne = length(find(F_it > 0));
fprintf(fid,'%d \n',ne) ; 
% Print out list of nodes for each
for k = 1:length(B)
    if F_it(k) > 0
        fprintf(fid,'%d \t %12.9f \n',k,F_it(k));
    end
end
fclose(fid);


