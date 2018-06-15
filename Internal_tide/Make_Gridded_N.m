%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Title:       Make_Gridded_N                                    %
%  Description: Read NETCDF of CTD casts with vertical profiles   %
%               of salinity temperature. Use TEOS toolbox to      %
%               compute profiles of the buoyancy frequency, N.    %
%  Inputs:      NETCDF files of gridded salinity and temeprature  %
%               data e.g. from the World Ocean Database:          %
%             https://www.nodc.noaa.gov/OC5/woa13/woa13data.html  %
%  Outputs:     A .mat file of profiles of N at each desired      %
%               contour level as specified at various locations   %
%  Requires:    TEOS toolbox http://www.teos-10.org/software.htm  %
%  Project:     Indian Ocean and Marginal Seas                    %
%  Author:      William Pringle                                   %
%  Created:     Oct 5 2016 , Updated Jun 15 2018                  %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;
close all; 

%% Set parameters here
sal = 'woa13_decav_s00_04v2.nc'; % salinity NETCDF file
temp  = 'woa13_decav_t00_04v2.nc'; % temperature NETCDF file

% Getting the 
z = ncread(sal,'depth');
lon = ncread(sal,'lon');
lat = ncread(sal,'lat');
SP = ncread(sal,'s_an');
T = ncread(temp,'t_an');
% Repeating z in longitude direction
[Z,~] = meshgrid(z,lon);

% Output filename
o_name      = 'Gridded_N_values'; %set .mat output filename (without suffix)

% Loop over data, convert Sp and T to Sa and CT, divide into bins,
% compute N2. 
L = length(lat);
z = 0.5*(z(1:end-1)+z(2:end));
N2 = zeros(length(lon),length(lat),length(z));
for i = 1:L
    if lat(i) > -60
    end
    SP_n = squeeze(SP(:,i,:))';
    T_n = squeeze(T(:,i,:))';
      
    % Convert the depth profile to a pressure profile
    p = gsw_p_from_z(-Z,lat(i));
    
    % Calculate absolute salinity from SP and p
    [SA, in_ocean] = gsw_SA_from_SP(SP_n,p',lon,lat(i));
    
    % Get the conservative temperature
    CT = gsw_CT_from_t(SA,T_n,p');
    
    % Calculate the profiles of N squared at p_mid
    [N2_n, p_mid] = gsw_Nsquared(SA,CT,p',lat(i));
    
    % Get the zbin levels back from p_mid
    zmid = gsw_z_from_p(p_mid,lat(i));
    
    N2(:,i,:) = N2_n';
end
N = sqrt(N2);
N = real(N);

%% Output results to a .mat file
save([o_name '.mat'],'lon','lat','z','N');
   