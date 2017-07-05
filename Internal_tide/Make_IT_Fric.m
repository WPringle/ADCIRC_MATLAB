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
type = 'directional';

% Use Nm or Nb for critical factor calculation
crit = 'Nb';

% Choose projection type (can be any, not restricted to evenly-gridded data)
proj = 'Mercator';
% Radius of earth for conversion to actual distances
R = 6378206.4; %[m]

% Load save data or not (set zero for first run for current grid, 
% set to 1 after running once for current grid)
load_data = 1;

% Coriolis coefficient
psi = 2*7.29212d-5;

% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency

% Choose the dimensionless coefficient for the internal tide friction
% (Only used if not Nycander tensor form)
C_it = 1.5;

% % Minimum and maximum allowable IT_Fric 
% MinFoldingTime = 6;   % hours
% MaxFoldingTime = 30;  % days
MinDepth         = 500; % m 

% ADCIRC mesh filame (fort.14)
adcirc_filename = '../INDWPAC_v1.8.mat';

% n data filename (only .mat file atm)
N_filename = 'E:\Global_Data\WOD_CTD_Casts\Global_N_values.mat';

% % original bathymetry filename if wanted for J calc
% bathyfile = ['E:\Global_Data\SRTM30_PLUS_w_Abyssal_Hills\' ...
%              'bathy_SSG_1_120_GLOBAL_landmask.nc'];

% load('E:\Global_Data\Nycander_IT_Fric\INDPAC_dJ_dh.mat',...
%      'dJx_g','dJy_g','lon_g','lat_g');
% dx = 1/120;
% lon = double(min(lon_g(:)):dx:max(lon_g(:)));
% lat = double(min(lat_g(:)):dx:max(lat_g(:)));
% [lon_g,lat_g] = ndgrid(lon,lat);
% Fx = griddedInterpolant(lon_g,lat_g,dJx_g','nearest');
% Fy = griddedInterpolant(lon_g,lat_g,dJy_g','nearest');
% dJx = Fx(VX); dJy = Fy(VX);
 
% save info filename
S_filename = 'fort14_it_info_direc';

%%%%%%%%%%Calculation region (should not need to alter below)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(adcirc_filename(end-2:end),'mat')
    load(adcirc_filename);
else
    [EToV,VX,B,~,~,~] = readfort14(adcirc_filename);
end
if load_data == 1
    %% Just load the data saved from previous run
    if strcmp(type,'scalar')
        load([S_filename '_directional.mat']);
    else
        load([S_filename '_' type '.mat']);
    end
else
    %% Load ADCIRC mesh and find the slope
    f = psi*sind(VX(:,2));           
    
    %% Load the constant contours of N values and compute Nb and Nmean
    load(N_filename); 
    
    [Nb,Nm] = Compute_Nb_Nm(VX(:,1),VX(:,2),B,zcontour,...
                            N,lon,lat,proj);  
                        
    %% Getting the J stuff if required (tensor type)
    if strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
        % Compute gradients of J from Nb, Nm and bathymetry B, and grid
        %[J,dJ,dh] = Compute_J_Nycander(EToV,VX,B,Nm,omega,...
        %                               5,MinDepth,proj,bathyfile,4);   
        [J,dJ,dh] = Compute_J_Nycander(EToV,VX,B,Nm,omega,...
                                       2,MinDepth,proj,[],4);  
        % save info for future computations 
        save([S_filename '_' type '.mat'],'VX','EToV','Nb','Nm','dh','f','B','J','dJ');                    
    else
        m_proj(proj,'lon',[ min(VX(:,1)) max(VX(:,1)) ],...
                    'lat',[ min(VX(:,2)) max(VX(:,2))]) 
        [xx,yy] = m_ll2xy(VX(:,1),VX(:,2));
        % Calculate slope directly on ADCIRC mesh
        [Hx,Hy] = ADCIRC_Bath_Slope( EToV,R*xx,R*yy, B );
        % save info for future computations 
        save([S_filename '_' type '.mat'],'Nb','Nm','Hx','Hy','f'); 
    end          
end
%
if MinDepth < 50
   I = knnsearch(VX(B >= 50,:),VX(B < 50 & B >= MinDepth,:)); 
   Nbn = Nb(B >= 50); Nmn = Nm(B >= 50); 
   Nb(B < 50 & B >= MinDepth) = Nbn(I); 
   Nm(B < 50 & B >= MinDepth) = Nmn(I); 
end

%if strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
%    H2_mesh = dh(:,1).^2 + dh(:,2).^2;      
%else
    H2_mesh = Hx.^2 + Hy.^2;
%end

%% Calculate F_it from Nb, Nm, Jx, Jy, and H2_mesh or as reqd.
if strcmp(type,'scalar')
   F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2)).*H2_mesh/omega;
elseif strcmp(type,'directional')
   F_it = C_it * sqrt((Nb.^2 - omega^2).*(Nm.^2 - omega^2))/omega;
elseif strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
   F_it = C_it * Nb/(4*pi).*sqrt(1-f.^2/omega^2)./B;
end
F_it = real(F_it); % in case becomes complex due to minus square root

%% Compute criticality and normalise the friction
if strcmp(crit(1:2),'Nm')
    alpha2  = (omega^2 - f.^2)./(Nm.^2 - omega^2);
elseif strcmp(crit(1:2),'Nb')
    alpha2  = (omega^2 - f.^2)./(Nb.^2 - omega^2);
end
if strcmp(crit(end),'f')
    alpha2 = str2double(crit(3:end-1))*alpha2;
end
if ~strcmp(crit,'none')
    gamma2 = max(H2_mesh./alpha2,1);
    % Normalise F_it by criticality
    if strcmp(crit(end),'k')
        % Use knife edge factor method
        % First divide by gamma
        %F_it = F_it./sqrt(gamma2);
        % where gamma > 1, multiply by 2 (knife edge to witch factor)
        %F_it(gamma2 > 1) = F_it(gamma2 > 1)*2;
        F_it = min(str2double(crit(end-1)),gamma2).*F_it./gamma2;
    else
        %gamma2(gamma2 > 1) = 2*gamma2(gamma2 > 1);
        F_it = F_it./gamma2;
    end
    % Make sure that if alpha2 < 0 that F_it is 
    % set equal to 0 since real part would be 0.
    F_it(alpha2 < 0) = 0;
end
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
% if strcmp(type,'tensor_to_scalar')
%    % Get the scalar value in terms of the directions of the tidal ellipse
%    F_it = F_it.*(dJ(:,1).*dh(:,1) + dJ(:,2).*dh(:,2));
%    L = length(find(F_it < 0));
%    F_it(F_it < 0) = 0;
% elseif strcmp(type,'tensor')
%     F_it(J == 0 | dh(:,1) == 0 | dh(:,2) == 0) = 0;
%     if ~strcmp(crit,'none')
%     for k = 1:length(F_it)
%         if F_it(k) > 0
%             % Output the C_11, C_12 = C_21, C_22 for the tensor
%             C_11 = 2*F_it(k)*dJ(k,1)*dh(k,1);
%             C_22 = 2*F_it(k)*dJ(k,2)*dh(k,2);
%             C_21 = F_it(k)*(dJ(k,1)*dh(k,2) + dJ(k,2)*dh(k,1));
%             A = [C_11 C_21; C_21 C_22];
%             if abs(min(eig(A))) > max(eig(A))
%                 % Ensure that the positive eigenvalue greater than negative
%                 % one
%                 %if any(abs(eig(A))/norm(A) > 1 + 1d-12)
%                 F_it(k) = 0;
%             end
%         end
%     end
%     end
% end
%
% In Antarctica lets make F_it greater since it is so small 
% (because of coarse resolution/bathy???)
%F_it(VX(:,2) < -64) = 4*4*F_it(VX(:,2) < -64);

%% Write out the fort.13 file
%if strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
%    filename = ['fort.13.' type '_' crit];
%else
    filename = ['fort.13.' type '_C' num2str(C_it) ...
                '_D' num2str(MinDepth) '_' crit] ;% '_Ant16'];
%end
fid = fopen(filename,'w');
% Print header
fprintf(fid,'%s \n','Spatial attributes description') ;
fprintf(fid,'%d \n',length(F_it)) ;   
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%s \n','internal_tide_friction') ;
fprintf(fid,'%s \n','1/time') ;
if strcmp(type,'tensor') || strcmp(type,'directional')
    fprintf(fid,'%d \n',3) ;
    fprintf(fid,'%f %f %f\n',0.0,0.0,0.0) ;  
else
    fprintf(fid,'%d \n',1) ;  
    fprintf(fid,'%f \n',0.0) ;  
end
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
            C_11 = 2*F_it(k)*dJx(k)*Hy(k);
            C_22 = 2*F_it(k)*dJy(k)*Hx(k);
            C_21 = F_it(k)*(dJx(k)*Hy(k) + dJy(k)*Hx(k));
            %C_11 = 2*F_it(k)*dJ(k,1)*dh(k,1);
            %C_22 = 2*F_it(k)*dJ(k,2)*dh(k,2);
            %C_21 = F_it(k)*(dJ(k,1)*dh(k,2) + dJ(k,2)*dh(k,1));
            fprintf(fid,'%d \t %15.9e %15.9e %15.9e \n',...
                    k,C_11,C_22,C_21);
        elseif strcmp(type,'directional')
            % Output the C_11, C_12 = C_21, C_22 for the tensor
            C_11 = F_it(k)*Hx(k)^2;
            C_22 = F_it(k)*Hy(k)^2;
            C_21 = F_it(k)*Hx(k)*Hy(k);
            fprintf(fid,'%d \t %15.9e %15.9e %15.9e \n',...
                    k,C_11,C_22,C_21); 
        else
            % Only need the F_it component
            fprintf(fid,'%d \t %15.9e \n',k,F_it(k));
            %fprintf(fid,'%d \t %15.9e \n',k,Nb(k));
        end
    end
end
fclose(fid);


