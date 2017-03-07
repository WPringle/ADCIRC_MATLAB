%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Title:       Make_IT_Fric                                      %
%  Description: Read constant contours of N values over depth     %
%               and bathymetry data to compute internal tide      %
%               friction coefficients on unstructured mesh        %
%  Inputs:      1) .mat files of N values at constant contours    %
%               2) Unstructured grid mesh with bathymetry         %                           
%  Outputs:     A fort.13 formatted file for use in ADCIRC/SMS    %
%  Project:     Indian Ocean and Marginal Seas                    %
%  Author:      William Pringle                                   %
%  Created:     Oct 5 2016                                        %
%  Updated:     Oct 24 2016, Dec 14 2016, Feb 25 2017             %
%  Requires:    functions - readfort14, Compute_Nb_Nm, m_proj,    %
%               m_ll2xy, ADCIRC_Bath_Slope, Compute_J_Nycander    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clearvars; 
clc; 
close all;
% 
%% Set parameters and filenames here
% Choose tensor (Nycander), directional or scalar form
type = 'tensor';

% Use Nm or Nb for critical factor calculation
crit = 'Nb';

% Choose projection type (can be any, not restricted to evenly-gridded data)
proj = 'Mercator';
% Radius of earth for conversion to actual distances
R = 6378206.4; %[m]

% Load save data or not (set zero for first run for current grid, 
% set to 1 after running once for current grid)
load_data = 0;

% Coriolis coefficient
psi = 2*7.29212d-5;

% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency

% Choose the dimensionless coefficient for the internal tide friction
% (Only used if not Nycander tensor form)
C_it = 1.5;

% % Minimum and maximum allowable IT_Fric 
% MinFoldingTime = 6;  % hours
% MaxFoldingTime = 30; % days
MinDepth       = 100; % m 

% ADCIRC mesh filame (fort.14)
adcirc_filename = '../IDIOMS_v7.1_SSG+TPXAnt_D2G.grd';

% n data filename (only .mat file atm)
N_filename = 'E:\Global_Data\WOD_CTD_Casts\Indian_Ocean_N_100m_processed.mat';

% original bathymetry filename if wanted for J calc
bathyfile = ['E:\Global_Data\SRTM30_PLUS_w_Abyssal_Hills\' ...
             'bathy_SSG_1_120_GLOBAL_landmask.nc'];

% save info filename
S_filename = 'fort14_it_info_SSG+TPXAnt.mat';

%%%%%%%%%%Calculation region (should not need to alter below)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    m_proj(proj,'lon',[ min(VX(:,1)) max(VX(:,1))],...
                'lat',[ min(VX(:,2)) max(VX(:,2))])
    [xx,yy] = m_ll2xy(VX(:,1),VX(:,2));             
    
    %% Load the constant contours of N values and compute Nb and Nmean
    load(N_filename); 
    %
    [Nb,Nm] = Compute_Nb_Nm(VX(:,1),VX(:,2),B,zcontour,...
                            N,lon,lat,proj);
    
%    load(S_filename);
                        
    %% Getting the J stuff if required (tensor type)
    if strcmp(type,'tensor')
        % Compute gradients of J from Nb, Nm and bathymetry B, and grid
        [J,dJ,dh] = Compute_J_Nycander(EToV,VX,B,Nm,omega,...
                                       5,MinDepth,proj,bathyfile,4);
        H2_mesh = dh(:,1).^2 + dh(:,2).^2;                          
        % save info for future computations 
        save(S_filename,'VX','EToV','Nb','Nm','dh','f','B','J','dJ');                    
    else
         % Calculate slope directly on ADCIRC mesh
        [Hx,Hy] = ADCIRC_Bath_Slope( EToV,R*xx,R*yy, B );
        H2_mesh = Hx.^2 + Hy.^2;
        % save info for future computations 
        save(S_filename,'VX','EToV','Nb','Nm','Hx','Hy','f','B'); 
    end          
end

%% Calculate F_it from Nb, Nm, Jx, Jy, and H2_mesh or as reqd.
if strcmp(type,'scalar')
   F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2)).*H2_mesh/omega;
elseif strcmp(type,'directional')
   F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2))/omega;
elseif strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
   F_it = Nb/(4*pi)*sqrt(1-f.^2/omega^2);
end
F_it = real(F_it); % in case becomes complex due to minus square root

%% Compute criticality and normalise the friction
if strcmp(crit,'Nm')
    alpha2  = (omega^2 - f.^2)./(Nm.^2 - omega^2);
elseif strcmp(crit,'Nb')
    alpha2  = (omega^2 - f.^2)./(Nb.^2 - omega^2);
end
gamma2 = max(H2_mesh./alpha2,1);

% Normalise F_it by criticality
F_it = F_it./gamma2;
% Make sure that if alpha2 < 0 that F_it is 
% set equal to 0 since real part would be 0.
F_it(alpha2 < 0) = 0;
%
% Clipping F_it for small and large e-folding times, and small depths
% MaxIT = 1/(MinFoldingTime*3600);
% MinIT = 1/(MaxFoldingTime*24*3600);
% 
% if strcmp(type,'scalar')
%     F_it(F_it < MinIT) = 0;
%     F_it(F_it > MaxIT) = MaxIT;
% elseif strcmp(type,'directional')
%     F_it(F_it.*H2_mesh < MinIT) = 0;
%     F_it(F_it.*H2_mesh > MaxIT) = MaxIT;
% end
F_it(B < MinDepth) = 0;
%
if strcmp(type,'tensor_to_scalar')
   % Get the scalar value in terms of the directions of the tidal ellipse
   F_it = F_it.*(dJ(:,1).*dh(:,1) + dJ(:,2).*dh(:,2));
   L = length(find(F_it < 0));
   F_it(F_it < 0) = 0;
end
%
%% Write out the fort.13 file
if strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
    filename = ['fort.13.' type '_' crit];
else
    filename = ['fort.13.' type '_C' num2str(C_it) '_' crit];
end
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
% Number of nodes not default (0.0)
ne = length(find(F_it > 0));
fprintf(fid,'%d \n',ne) ; 
% Print out list of nodes for each
for k = 1:length(F_it)
    if F_it(k) > 0
        if strcmp(type,'tensor')
            % Output the C_11, C_12 = C_21, C_22 for the tensor
            C_11 = 2*F_it(k)*dJ(k,1)*dh(k,1);
            C_22 = 2*F_it(k)*dJ(k,2)*dh(k,2);
            C_21 = F_it(k)*(dJ(k,1)*dh(k,2) + dJ(k,2)*dh(k,1));
            fprintf(fid,'%d \t %15.9e %15.9e %15.9e \n',...
                    k,C_11,C_22,C_21);
        else
            % Only need the F_it component
            fprintf(fid,'%d \t %15.9e \n',k,F_it(k));
        end
    end
end
fclose(fid);


