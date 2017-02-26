%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Title:       Make_N_Profile                                    %
%  Description: Read NETCDF of CTD casts with vertical profiles   %
%               of salinity temperature. Use TEOS toolbox to      %
%               compute profiles of the buoyancy frequency, N.    %
%  Inputs:      NETCDF files of CTD casts from the World Ocean    %
%               Database. Obtain from:                            %
%     http://www.nodc.noaa.gov/OC5/SELECT/dbsearch/dbsearch.html  %
%  Outputs:     A .mat file of profiles of N at each desired      %
%               contour level as specified at various locations   %
%  Requires:    TEOS toolbox http://www.teos-10.org/software.htm  %
%  Project:     Indian Ocean and Marginal Seas                    %
%  Author:      William Pringle                                   %
%  Created:     Oct 5 2016                                        %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc;
close all; 

%% Set parameters here
filename = 'ocldb1475086690.12161.CTD.nc'; % filename of master NETCDF file
%ncdisp(filename) % use to display contents of NETCDF if required
cast_direc = 'Casts/'; % directory of all the casts
cast_prefix = 'wod_';  % These should be as specified here once downloaded
cast_suffix = 'O.nc';  % and extracted

% Define your region of interest
lat_min  = 10; %Set min latitude of bounding box
lat_max  = 20; %Set max latitude of bounding box
lon_min  = 80; %Set min longitude of bounding box
lon_max  = 90; %Set max longitude of bounding box

% Define the bin intervals to get N profiles at
% Note that bin intervals too small will give very noisy N profiles
bin_int = 100; %Default is 100 m as suggested by:
% King, B., et al. (2012). Buoyancy frequency profiles and internal 
% semidiurnal tide turning depths in the oceans. JGR: Oceans, 117(C4)
% 100 m is useful for looking at N over great depths well below
% the thermocline. For looking inside thermocline in detail this should be 
% set much smaller, e.g. 10 m (and maximum depth contour set smaller too)

% Maximum depth contour to output to (contours: bin_int:bin_int:max_depth)
max_depth = 1d4; 

% Quality control parameters
Point_length_min = 10; % Minimum number of data points per bin
num_std          = 3;  % Number of standard deviations allowable from 
                       % spatial mean of neighbouring points
search_radius    = 5;  % degree search radius for elimination above

% Output filename
o_name      = 'SCS_N_values'; %set .mat output filename (without suffix)

% Set if you want to plot some casts to check goodness of calc
plot_to_check = 1;  %0 - no plot, 1 - plot 
plot_num      = 10; % number of plots to see if plot_to_check = 1

% Set if you want to plot some scatter points
plot_scat     = 1;  %0 - no plot, 1 - plot  
plot_num1     = 7;  % number of plots to see if plot_scat = 1

%%
%%%%%%%%%%Calculation region (should not need to alter below)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read in the cast IDs, their locations and time of observation
cast = ncread(filename,'cast');
lat  = ncread(filename,'lat');
lon  = ncread(filename,'lon');
time = ncread(filename,'time');
% 
% Find bounds of the region
I = find(lat < lat_min | lat > lat_max | ...
         lon < lon_min | lon > lon_max);
% Elminate data outside bounds
lat(I) = []; lon(I) = []; cast(I) = []; time(I) = [];
% 
% Extract the profiles of z, T and SP
z  = cell(length(cast),1); T = cell(length(cast),1); 
SP = cell(length(cast),1);
for c = 1:length(cast)
    cc = num2str(cast(c),'%09d');
    cast_n = [cast_direc cast_prefix cc cast_suffix];
    
    ncid = netcdf.open(cast_n,'nowrite'); flag = 0;
    % Try Salinity & Temperature variables to see if exist
    for VarName = {'Salinity','Temperature'}
        try
            ID = netcdf.inqVarID(ncid,char(VarName));
        catch exception
            if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                flag = 1;
                break;
            end
        end
    end
    netcdf.close(ncid)
    % Skip if salinity or temperature does not exist
	if flag == 1; continue; end; 

    % Get the z profiles
    z{c} = ncread(cast_n,'z');
    % Get the temperature profiles
    T{c} = ncread(cast_n,'Temperature');
    % Get the practical salinity profiles
    SP{c} = ncread(cast_n,'Salinity');
end
% Eliminate the empty data 
%# find empty cells
emptyCells = cellfun('isempty',SP);
%# remove empty cells
SP(emptyCells) = []; z(emptyCells) = []; T(emptyCells) = [];
lat(emptyCells) = []; lon(emptyCells) = []; time(emptyCells) = [];
% 
%
%% Initialise the N2 (N squared) and zbins
N2 = cell(size(lat)); zbins = cell(size(lat));
% Loop over data, convert Sp and T to Sa and CT, divide into bins,
% compute N2. 
L = length(lat);
for i = 1:L;
    % Skip if quality sucks
    % Only get data with good amount of points per bin and 
    if z{i}(end) < bin_int; continue; end
    if length(z{i}) < z{i}(end)*Point_length_min/bin_int; continue; end
    % Eliminate points deeper than usually physically possible (bad data?)
    if z{i}(end) > 1d4;
        I = find(z{i} > 1d4);
        SP{i}(I) = []; T{i}(I) = []; z{i}(I) = []; 
    end
    % Convert the depth profile to a pressure profile
    p = gsw_p_from_z(-z{i},lat(i));
    
    if any(p < -1.5 | p > 12000)
        continue;
    end
    
    % Calculate absolute salinity from SP and p
    [SA, in_ocean] = gsw_SA_from_SP(SP{i},p,lon(i),lat(i));
    
    % Skip if your point is not in the ocean
    if in_ocean(1) == 0; continue; end
    
    % Skip if quality of SA sucks
    if length(SA(SA > 0)) < z{i}(end)*Point_length_min/bin_int; continue; end
    
    % Get the conservative temperature
    CT = gsw_CT_from_t(SA,T{i},p);
    
    % Make bins every bin_int metres 
    binz = 0:bin_int:z{i}(end);   
    
    % Convert z bins to pressure
    binp = gsw_p_from_z(-binz,lat(i));
    % Make sure last p point is included
    binp(end+1) = p(end);
    % Process the pressure into the pressure bins
    [B,idx] = histc(p,binp);
    
    % Evaluate CT and SA at the pressure bins
    CTb = accumarray(idx(idx > 0),CT(idx > 0),[],@mean);
    SAb = accumarray(idx(idx > 0),SA(idx > 0),[],@mean);
    
    % Get the middle points of the bins
    binpm = 0.5*(binp(1:end-1) + binp(2:end));
    
    % Calculate the profiles of N squared at p_mid
    [N2{i}, p_mid] = gsw_Nsquared(SAb(1:end-1),CTb(1:end-1),binpm',lat(i));
    
    % Get the zbin levels back from p_mid (should be intervals of bin_int
    % except when z{i}(end) is < 2*bin_int
    zbins{i} = gsw_z_from_p(p_mid,lat(i));
end
% Eliminate the empty data 
%# find empty cells
emptyCells = cellfun('isempty',N2);
%# remove empty cells
N2(emptyCells) = []; zbins(emptyCells) = [];
lat(emptyCells) = []; lon(emptyCells) = []; time(emptyCells) = [];

% Check data if required by plotting each n
L = length(lat);
for i = 1:L
    if plot_to_check == 1
        if mod(i,floor(L/plot_num)) == 0 
            figure;
            plot(N2{i},zbins{i},'o-')
        end 
    end
end

%% Process the data into desired contours eliminating bad quality data
% Make the contours
zcontour = bin_int:bin_int:max_depth;
% Put lat lon into latN lonN
latN = double(lat);
lonN = double(lon);

% Initialise the N and output lat and lon
N = cell(size(zcontour));
lon = cell(size(zcontour));
lat = cell(size(zcontour));
for zvalue = zcontour
    
    % Eliminate values where N2 does not exist at current contour level
    N2z = zeros(length(N2),1);
    for i = 1:length(N2)
        if length(N2{i}) >= zvalue/bin_int
            N2z(i) = N2{i}(zvalue/bin_int);     
        end
    end
    lonNz = lonN(N2z > 0);
    latNz = latN(N2z > 0);
    N2z  = N2z(N2z > 0);

    % Eliminate shitty values outside of 
    % specified standard deviations for specified degree radius search
    IDX = rangesearch([lonNz, latNz],[lonNz, latNz],search_radius);
    for i = 1:length(lonNz)
       N2mean = mean(N2z(IDX{i}));
       N2std  = std(N2z(IDX{i}));
       N2var = abs(N2z(i) - N2mean);
       if N2var > num_std*N2std;
            N2z(i) = 0;
       end
    end
    lonNz = lonNz(N2z > 0);
    latNz = latNz(N2z > 0);
    N2z  = N2z(N2z > 0);

    % Make the N, lon, and lat profiles
    lon{zvalue/bin_int} = lonNz;
    lat{zvalue/bin_int} = latNz;
    N{zvalue/bin_int}   = sqrt(N2z); 
end   
%% Output results to a .mat file
save([o_name '.mat'],'lon','lat','N','zcontour');
   
%% plot scatter at contours if required
if plot_scat == 1;
    zc = linspace(zcontour(1),zcontour(end),plot_num1);
    zc = bin_int*int32(zc/bin_int);
    for zvalue = zc
        i_now = zvalue/bin_int;
        if ~isempty(N{i_now})
            figure;
            scatter(lon{i_now},lat{i_now},[],N{i_now},'filled')
            title(['z =' num2str(zvalue) ' m']) 
            colorbar
            caxis([0 prctile(N{i_now},99)])
        end
    end
end