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
%  Updated:     Oct 24 2016, Dec 14 2016                          %
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
lon0 = 75.015117 * pi/180; lat0 = -30.890245 * pi/180;
R = 6378206.4;

% Choose tensor or scalar form
type = 'tensor';

% Use Nm or Nb for critical factor calculation
crit = 'Nm';

% Load save data or not (set zero for first run for current grid, 
% set to 1 after running once for current grid)
load_data = 0;

% Coriolis coefficient
psi = 2*7.29212d-5;

% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency

% Choose the dimensionless coefficient for the internal tide friction
C_it = 1/2;

% Minimum and maximum allowable IT_Fric 
MinFoldingTime = 6;  % hours
MaxFoldingTime = 30; % days
MinDepth       = 100; % m 

MaxIT = 1/(MinFoldingTime*3600);
MinIT = 1/(MaxFoldingTime*24*3600);

% ADCIRC mesh filame (fort.14)
adcirc_filename = 'fort.14';

% n data filename (only .mat file atm)
N_filename = 'E:\Global_Data\WOD_CTD_Casts\Indian_Ocean_N_100m_processed.mat';

% save info filename
S_filename = 'fort14_it_info.mat';
%%
%%%%%%%%%%Calculation region (should not need to alter below)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deprecated section calculating slope directly from GEBCO 
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

if load_data == 1
    %% Just load the data saved from previous run
    load(S_filename);
    H2_mesh = Hx.^2 + Hy.^2;
else
    %% Load ADCIRC mesh and find the slope
    if strcmp(adcirc_filename(end-2:end),'mat')
        load(adcirc_filename);
    else
        [EToV,VX,B,~,~,~] = readfort14(adcirc_filename);
    end
    f = psi*sind(VX(:,2));
    
    xx = VX(:,1) * pi/180;  yy = VX(:,2)* pi/180; 
    xx = R * (xx - lon0) * cos(lat0);
    yy = R * yy;

    % Do the interpolation onto the ADCIRC mesh
    [Hx,Hy] = ADCIRC_Bath_Slope( EToV,xx,yy, B );
    H2_mesh = Hx.^2 + Hy.^2;
    %
    %% Load the constant contours of N values and compute Nb and Nmean
    load(N_filename); 
    %
    [Nb,Nm] = Compute_Nb_Nm(VX(:,1),VX(:,2),B,zcontour,...
                            N,lon,lat,lon0,lat0);
    
    % save info for future computations 
    save(S_filename,'Nb','Nm','Hx','Hy','f');                  
end

%% Calculate F_it from Nb, Nm, and H2_mesh
if strcmp(type,'scalar')
   F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2)).*H2_mesh/omega;
elseif strcmp(type,'tensor')
   F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2))/omega;
end
F_it = real(F_it); % in case becomes complex due to minus square root

%% Compute criticality and normalise the friction
if strcmp(crit,'Nm')
    alpha2  = (omega^2 - f.^2)./(Nm.^2 - omega^2);
elseif strcmp(crit,'Nn')
    alpha2  = (omega^2 - f.^2)./(Nn.^2 - omega^2);
end
gamma2 = max(H2_mesh./alpha2,1);

% Normalise F_it by criticality
F_it = F_it./gamma2;
% Make sure that if alpha2 < 0 that F_it is 
% set equal to 0 since real part would be 0.
F_it(alpha2 < 0) = 0;
% 
%% Write out the fort.13 file
filename = ['fort.13.' type '_C' num2str(C_it) '_' crit];
fid = fopen(filename,'w');
% Print header
fprintf(fid,'%s \n','Spatial attributes description') ;
fprintf(fid,'%d \n',length(F_it)) ;   
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%s \n','internal_tide_friction') ;
fprintf(fid,'%s \n','1/time') ;
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%f \n',0.0) ;  
fprintf(fid,'%s \n','internal_tide_friction') ;
% 
% Clipping F_it for small and large e-folding times, and small depths
%F_it(F_it > MaxIT) = MaxIT;
if strcmp(type,'scalar')
    F_it(F_it < MinIT) = 0;
elseif strcmp(type,'tensor')
    F_it(F_it.*H2_mesh < MinIT) = 0;
end
%F_it(B < MinDepth) = 0;
%
% Number of nodes not default (0.0)
ne = length(find(F_it > 0));
fprintf(fid,'%d \n',ne) ; 
% Print out list of nodes for each
for k = 1:length(F_it)
    if F_it(k) > 0
        fprintf(fid,'%d \t %12.9f \n',k,F_it(k));
    end
end
fclose(fid);


