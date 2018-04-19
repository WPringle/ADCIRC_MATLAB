function [dJx, dJy, lon_c, lat_c] = ...
               Calc_dJ_Nyc_Struc(lon,lat,bathy,Nm,omega,hcut,acut,bbox,par)
% function [dJx, dJy, lon_c, lat_c] =   
%              Calc_dJ_Nyc_Struc(lon,lat,bathy,Nm,omega,hcut,acut,bbox,par) 
%
% Calc_dJ_Nyc_struc: Computes the gradients of J for use in the 
% Nycander (2005) wave drag formulation on a unstructured grid using the 
% Green and Nycander (2013) FD method
%
% Inputs    :   lon   - vector or meshgrid matrix of longitudes
%               lat   - vector or meshgrid matrix of latitudes
%               bathy - gridded bathymetric data (minus of the depth)
%               Nm    - gridded depth-averaged buoyancy frequency
%               omega - the circular frequency (rad/s) of the tidal 
%                       component of interest (usually the M2 one)
%               hcut  - The cutoff height (only find dJx.. when bathy < hcut)
%               acut  - The cutoff radius (only find dJx.. when a > acut)
%               bbox  - the bounding box of lon, lat to calculate dJ within
%                       Set bbox = [lon_min, lon_max; lat_min, lat_max]
%               par   - <= 1, do not perform parallel for loop
%                       > 1, parallel for loop where the par integer 
%                       indicates the number of processes to use
%
% Outputs    :  dJx   - the slope of J in x direction
%               dJy   - the slope of J in y direction 
%               lon_c - meshgrid of longitudinal points where dJ is defined
%               lat_c - meshgrid of latitudinal points where dJ is defined
%
% Author     : William Pringle, Dec 12 2017
%
% Requires   : 'm_map'toolbox for calculating distances between points 
%              found at https://www.eoas.ubc.ca/~rich/map.html
%
% References : Nycander, J. (2005). Generation of internal waves in the deep
%              ocean by tides. Journal of Geophysical Research, 
%              110(C10), C10028. http://doi.org/10.1029/2004JC002487
%
%              Green, J. A. M., & Nycander, J. (2013). A Comparison of Tidal
%              Conversion Parameterizations for Tidal Models. Journal of 
%              Physical Oceanography, 43(1), 104–119. 
%              http://doi.org/10.1175/JPO-D-12-023.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting constants
% Coriolis coefficient
psi = 2*7.29212d-5;
% Beta coefficient (Nycander, 2005) in radius
beta = 1.455;
% Cutoff factor in calculation of dJ
rcutfac = 5 ; % points within |r - r'| < rcutfac*a will be used
                
%% Setting up parallel if desired
if par > 0
    if isempty(gcp())
        parpool(par);
    end
end

%% Checking inputs
% Check if lon lat are vectors or meshgrid (make into meshgrid if latter)
if any(size(lon) == 1)
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
if size(Nm,1) ~= size(lon,1)
   Nm = Nm';
   if size(Nm,1) ~= size(lon,1)
       disp('Nm is badly sized')
       return;
   end
end

%% Evaluate bathy and buoyancy frequency at the midpoints of the input 
%% nodes to evaluate dJ at, and remove values outside of the bounding box
% Get midpoints coordinates
lat_c = 0.5*(lat(1:end-1,1:end-1) + lat(2:end,1:end-1));
lon_c = 0.5*(lon(1:end-1,1:end-1) + lon(1:end-1,2:end));

% Get the bathymetry and buoyancy frequencies at the midpoints 
% of the four surrounding nodes
h_c   = 0.5*(bathy(1:end-1,:) + bathy(2:end,:));
h_c   = 0.5*(h_c(:,1:end-1) + h_c(:,2:end));
Nm_c  = 0.5*(Nm(1:end-1,:) + Nm(2:end,:));
Nm_c  = 0.5*(Nm_c(:,1:end-1) + Nm_c(:,2:end));

% Find bounds
I = find( lon_c(1,:) > bbox(1,1) & lon_c(1,:) <= bbox(1,2));
J = find( lat_c(:,1) > bbox(2,1) & lat_c(:,1) <= bbox(2,2));

% Only keep part we want to calc.
lat_c = lat_c(J,I); lon_c = lon_c(J,I); h_c = h_c(J,I); Nm_c = Nm_c(J,I);

%% Computing the cutoff length of the Green's function
f = psi*sind(lat_c) ;            % Coriolis coefficients
dom = pi*sqrt(omega^2 - f.^2);   % dominator of the cutoff length
a = (beta./dom).*-h_c.*Nm_c ;    % the cutoff length
a = real(a);                     % just take real part because dom is 
                                 % imaginary above some latitude

% Get the d_x values on the wgs84 Spheroid
d_x = m_idist(lon(J(1):J(end),I(1):I(end)),lat(J(1):J(end),I(1):I(end)),...
            lon(J(1):J(end),I(2):I(end)+1),lat(J(1):J(end),I(2):I(end)+1));
          
% Get the d_y values on the wgs84 Spheroid
d_y = m_idist(lon(J(1):J(end),I(1):I(end)),lat(J(1):J(end),I(1):I(end)),...
            lon(J(2):J(end)+1,I(1):I(end)),lat(J(2):J(end)+1,I(1):I(end)));

% number of nodes to search back and forwards
y_n = ceil(rcutfac*a./d_y);
x_n = ceil(rcutfac*a./d_x);

% Initialising the arrays
dJx = zeros(size(lat_c)) ;
dJy = zeros(size(lat_c)) ;

% Compute the singlularity correction coefficients
dh = bathy(J+1,I+1) + bathy(J,I+1) - bathy(J+1,I) - bathy(J,I);
s = d_x./d_y;
correc = ((2./s).*log(sqrt(s.^2 +1) + s) - 4*s.^2./(1+s.^2).^(3/2)).*dh;

lat_n = length(J);
lon_n = length(I);

dhx = zeros(size(correc));
dhy = zeros(size(correc));

for l_y = 1:lat_n
    disp(['l_y = ' num2str(l_y)])
    % Get global ID no.s
    ly_g = J(l_y);
    for l_x = 1:lon_n
        
        % Skip for small depths and small radii
        if h_c(l_y,l_x) > hcut; continue; end
        if a(l_y,l_x) < acut; continue; end
        
        % Get global ID no.s
        lx_g = I(l_x);
        
        % Get the slopes at the center point
        dhx(l_y,l_x) = 0.5*(bathy(ly_g+1,lx_g+1) + bathy(ly_g,lx_g+1) ...
           - bathy(ly_g+1,lx_g) - bathy(ly_g,lx_g))/d_x(l_y,l_x);
        dhy(l_y,l_x) = 0.5*(bathy(ly_g+1,lx_g+1) + bathy(ly_g+1,lx_g) ...
           - bathy(ly_g,lx_g+1) - bathy(ly_g,lx_g))/d_y(l_y,l_x); 
        
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
        bathy_n = bathy(range_y,range_x);
        dh = bathy_n(K) - h_c(l_y,l_x);
        
        % Get the gradient of greens' function
        dg = diffgrnfunc( r(K), a(l_y,l_x) );
        
        % The non-directional component of the integral
        NDI = dh.*dg./r(K);
        
        % Compute the sums with the dx to get the integral
        Jx_sum = sum(NDI.*dx_r);
        Jy_sum = sum(NDI.*dy_r);

        % Multiply by dx,dy and add the correction for the singularity
        dJx(l_y,l_x) = Jx_sum * d_x(l_y,l_x)*d_y(l_y,l_x) + correc(l_y,l_x);
        dJy(l_y,l_x) = Jy_sum * d_x(l_y,l_x)*d_y(l_y,l_x) + correc(l_y,l_x);     
    end
    disp(['dJx = ' num2str(dJx(l_y,l_x))])
end
if par > 0
    delete(gcp('nocreate'));
end
% Transpose outputs if the inputs were the transpose of the meshgrid 
if Tr
   dJx = dJx'; dJy = dJy';
end
%EOF
end

% Some stuff could use to try speed up; not done
% Get the coordinate Matrices
%xx = cumsum(d_x,2); 
%yy = cumsum(d_y,1); 
% Find where we need to calc
%idv = find( h_c < hcut  & a > acut );
% % Make the coordinate matrix
% XX = [xx(idv) yy(idv)];
% % set up kdtree
% Mdl = KDTreeSearcher(XX); 
% % sort by radius
% [aiv, ia] = sort(a(idv));
% xiv = XX(ia,:);
% % number of times to do rangesearch
% nj = length(idv) ;
% blocks = ceil(nj*rcutfac^2*1.5e-6);
% nn = 0;
% disp(['number of blocks = ' num2str(blocks)])
% for ii = 1:blocks
%     tic
%     ns = int64((ii-1)*nj/blocks)+1;
%     ne = int64(ii*nj/blocks);
%     % Do range search
%     [Idx, Dist] = rangesearch(Mdl,xiv(ns:ne,:),rcutfac*aiv(ne));  
%     % Loop through Idx
%     for jj = 1:length(Idx)
%         nn = nn + 1;
%         val = gaussC(Dist{jj}(Dist{jj} < rcutfac*aiv(nn)), aiv(nn)^2); % Get the weights
%         val = val/norm(val,1); % divide by the norm
%         %val2 = greensF(Dist{jj}(Dist{jj} < rcutfac*aiv(nn))/aiv(nn)); % Get the weights
%         NbTemp(nn) = val*Nb(Idx{jj}(Dist{jj} < rcutfac*aiv(nn))); %sum the Nb with weights
%         % NbTemp(nn) = mean(Nb(Idx{jj}(Dist{jj} < aiv(nn))));
%     end
%     disp(['finished block no. = ' num2str(ii)])
%     toc
% end
% 
% mcell = max(max(x_n(:)));
% % loop over x-dir searching nodes
% for xx = 1:mcell
%     % Get the difference of the topographic heights
%     JJ      = find(x_n == xx);
%     [r,c]   = ind2sub(size(lon_c),JJ);
%     ly_g = J(c);
%     lx_g = I(r);
%     range_x = max(lx_g - xx+1,1):min(lx_g + xx,size(lon,2));
%     range_y = max(ly_g - y_n(JJ)+1,1):min(ly_g + y_n(JJ),size(lat,1));
%     
%     back = bathy(J,I-xx+1);
%     forw = bathy(J,I+xx);
%     dhb = back(JJ) - h_c(JJ); 
%     dhf = forw(JJ) - h_c(JJ); 
% end