function [dJx, dJy, dhx, dhy, lon_c, lat_c] = ...
            dJ_Nycander_grid(lon,lat,bathy,N_file,omega,hcut,bbox,parallel)
% Computes dJ for use in Nycander (2005) wave drag
% formula on the original bathymetric grid using FD fomulation
% function [dJx, dJy]= dJ_Nycander_grid(bathyfile,Nm,omega,hcut,proj,parallel)   
% 
% Inputs:  lon  - vector or meshgrid matrix of longitudes
%          lat  - vector or meshgrid matrix of latitudes
%         bathy - gridded bathymetric data (minus of the depth)
%        N_file - .mat file of the contour info of the Brunt–Väisälä 
%                 buoyancy frequency 
%         omega - the circular frequency (rad/s) of the tidal component of
%                 interest (usually the M2 one)
%          hcut - The cutoff height (only find dJx.. when bathy > hcut)
%          bbox - the bounding box of lon, lat to calculate dJx, dJy in.
%                 set; bbox = [lon_min, lon_max; lat_min, lat_max]
%      parallel - perform parallel for loop or not (>0, yes where the integer
%                 indicates number of processes to use, = 0, no)
%
% Outputs: dJx - the slope of J in x direction
%          dJy - the slope of J in y direction 
%          dhx - the slope of bathy in x direction
%          dhy - the slope of bathy in y direction 
%        lon_c - meshgrid of longitudinal points where dJ, dh are defined
%        lat_c - meshgrid of latitudinal points where dJ, dh are defined
%
% Authors : William Pringle March 13 2017
%
% Requires : - m_map toolbox for projection and calculating differences
%              between pointd

% Reference: For derivation and fomula of the wave drag - 
%            Nycander, J. (2005). Generation of internal waves in the deep
%            ocean by tides. Journal of Geophysical Research, 
%            110(C10), C10028. http://doi.org/10.1029/2004JC002487
%
%            For the "improved" finite-difference formula
%            Green, J. A. M., & Nycander, J. (2013). A Comparison of Tidal
%            Conversion Parameterizations for Tidal Models. Journal of 
%            Physical Oceanography, 43(1), 104–119. 
%            http://doi.org/10.1175/JPO-D-12-023.1

% Set some constants
% Coriolis coefficient
psi = 2*7.29212d-5;
% beta coefficient
beta = 1.455;

% Some arbitrary settings for the numerics 
rcutfac = 2.5 ; % Elements within |r - r'| < rcutfac*a will be used in the calculation of dJ/dx
acut = 2   ; % if a < acut (in km), ignore 

% Setting up parallel
if parallel > 0
    delete(gcp('nocreate'));
    parpool(parallel);
end

%% Checking inputs
% Check if lon lat are vectors or meshgrid (make into meshgrid if latter)
if any(size(lon)) == 1
    [lon, lat] = meshgrid(lon,lat);
end
% Check the size of bathy and Nm is correct
Tr = 0;
if size(bathy,1) ~= size(lon,1)
   bathy = bathy';
   Tr = 1;
   if size(bathy,1) ~= size(lon,1)
       disp('bathy is badly sized')
       return;
   end
end

% Get midpoints
lat_c = 0.5*(lat(1:end-1,1:end-1) + lat(2:end,1:end-1));
lon_c = 0.5*(lon(1:end-1,1:end-1) + lon(1:end-1,2:end));

% Get the midpoints of the four surrounding topographic points
h_c   = 0.5*(bathy(1:end-1,:) + bathy(2:end,:));
h_c   = 0.5*(h_c(:,1:end-1) + h_c(:,2:end));

% Find bounds
I = find( lon_c(1,:) > bbox(1,1) & lon_c(1,:) < bbox(1,2));
J = find( lat_c(:,1) > bbox(2,1) & lat_c(:,1) < bbox(2,2));

% Only keep part we want to calc.
lat_c = lat_c(J,I); lon_c = lon_c(J,I); h_c = h_c(J,I); 

%% Get Nm at the midpoints
Bf = load(N_file); 
%
[~,Nm] = Compute_Nb_Nm(lon_c,lat_c,-h_c,Bf.zcontour,...
                       Bf.N,Bf.lon,Bf.lat,'Mercator');

%% Computing the cutoff length of green's function
f = psi*sind(lat_c) ;            % Coriolis coefficients
dom = pi*sqrt(omega^2 - f.^2);   % dominator of the cutoff length
a = (beta./dom).*-h_c.*Nm ;      % the cutoff length

% Get the d_x values on the wgs84 spheroid
d_x = m_idist(lon(J(1):J(end),I(1):I(end)),lat(J(1):J(end),I(1):I(end)),...
              lon(J(1):J(end),I(2):I(end)+1),lat(J(1):J(end),I(2):I(end)+1));
% Get the d_y values on the wgs84 spheroid
d_y = m_idist(lon(J(1):J(end),I(1):I(end)),lat(J(1):J(end),I(1):I(end)),...
              lon(J(2):J(end)+1,I(1):I(end)),lat(J(2):J(end)+1,I(1):I(end)));

% number of nodes to search back and forwards
y_n = ceil(rcutfac*a./d_y);
x_n = ceil(rcutfac*a./d_x);

lat_n = size(lat_c,1);
lon_n = size(lon_c,2);

% Initialising the arrays
dJx = zeros(size(lat_c)) ;
dJy = zeros(size(lat_c)) ;
dhx = zeros(size(lat_c)) ;
dhy = zeros(size(lat_c)) ;
parfor l_y = 1:lat_n
    disp(['l_y = ' num2str(l_y)])
    % Get global ID no.s
    ly_g = J(l_y);
    for l_x = 1:lon_n
        
        % Skip for small depths and small radii
        if h_c(l_y,l_x) > hcut; continue; end
        if a(l_y,l_x)/1000 < acut; continue; end
        
        % Get global ID no.s
        lx_g = I(l_x);
        
        % Get the slopes at the center point
        dhx(l_y,l_x) = 0.5*(bathy(ly_g+1,lx_g+1) + bathy(ly_g,lx_g+1) ...
           - bathy(ly_g+1,lx_g) - bathy(ly_g,lx_g))/d_x(l_y,l_x);
        dhy(l_y,l_x) = 0.5*(bathy(ly_g+1,lx_g+1) + bathy(ly_g+1,lx_g) ...
           - bathy(ly_g,lx_g+1) - bathy(ly_g,lx_g))/d_y(l_y,l_x); 
       
        % Get the correction coefficient
        dh = bathy(ly_g+1,lx_g+1) + bathy(ly_g,lx_g+1) ...
           - bathy(ly_g+1,lx_g) - bathy(ly_g,lx_g);
        s = d_x(l_y,l_x)/d_y(l_y,l_x);
        correc = ((2/s)*log(sqrt(s^2 +1) + s) - 4*s^2/(1+s^2)^(3/2))*dh;
        
        % Get the range to search
        range_x = max(lx_g - x_n(l_y,l_x)+1,1):min(lx_g + x_n(l_y,l_x),size(lon,2));
        range_y = max(ly_g - y_n(l_y,l_x)+1,1):min(ly_g + y_n(l_y,l_x),size(lat,1));
        
        % Compute the distances
        [r,a12] = m_idist(lon_c(l_y,l_x),lat_c(l_y,l_x),...
                    lon(range_y,range_x),lat(range_y,range_x));
        K = find(r < rcutfac*a(l_y,l_x));              
        %
        if isempty(K)
            continue;
        end
        % Get x and y components of the distance
        dx_r = -sind(a12(K)).*r(K);
        dy_r = -cosd(a12(K)).*r(K);
        
        % Get the difference of the topographic heights
        bathy_n = bathy(range_x,range_y);
        dh = bathy_n(K) - h_c(l_y,l_x);
        
        % Get the gradient of greens' function
        dg = diffgrnfunc( r(K), a(l_y,l_x) );
        
        % The non-directional component of the integral
        NDI = dh.*dg./r(K);
        
        % Compute the sums with the dx to get the integral
        Jx_sum = sum(NDI.*dx_r);
        Jy_sum = sum(NDI.*dy_r);

        % Multiply by dx,dy and add the correction for the singularity
        dJx(l_y,l_x) = Jx_sum * d_x(l_y,l_x)*d_y(l_y,l_x) + correc;
        dJy(l_y,l_x) = Jy_sum * d_x(l_y,l_x)*d_y(l_y,l_x) + correc;
    end
end
if parallel > 0
    delete(gcp('nocreate'));
end
%EOF
end
