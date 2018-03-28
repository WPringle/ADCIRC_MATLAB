function [J,dJ] = Calc_dJ_Nyc_Unstruc(EToV,VX,B,Nm,omega,...
                                      pcub_c,hcut,proj,parallel)
% Computes J, and its derivatives for use in Nycander (2005) wave drag
% formula on the finite-element grid using quadrature of arbitrary degree
% function [J,dJ] = Compute_J_Nycander(EToV,VX,B,Nm,omega,pcub_c,hcut,proj,parallel)
% Inputs:  EToV - the element table (Num_ele x 3)     
%          VX   - the vertex positions (Num_node x 2)
%          B    - the bathymetric depth at each node (Num_node x 1)
%          Nm   - the depth-averaged Brunt–Väisälä buoyancy frequency at
%                 each node (Num_node x 1)
%         omega - the circular frequency (rad/s) of the tidal component of
%                 interest (usually the M2 one)
%        pcub_c - the degree of the quadrature required (set = 2 if unsure)
%          hcut - The cutoff depth (only find J.. when H > hcut)
%          proj - a string defining the projection to use, e.g. 'Mercator'.
%                 See m_proj('get') for other options in m_map 
%      parallel - perform parallel for loop or not (>0, yes where the integer
%                 indicates number of processes to use, = 0, no)
%
% Outputs: J  - The J variable (Num_node x 1)
%          dJ - The slope of J in both directions (dJ/dx and dJ/dy)
%               (Num_node x 2)
%
% Authors : Damrongsak Wirasaet, Feb 2017 - did all calculation stuff
%           William Pringle Feb 27 2017, just made into function to
%           interface with Make_IT_Fric
%
% Requires : - Suite of functions for quadrature calc. found in
%             "Internal_tide/Quad_functions" folder on ADCIRC_MATLAB github
%            - m_map toolbox for projection
%            (https://www.eoas.ubc.ca/~rich/map.html)
%            - Distmesh toolbox for dcircle function
%            (http://persson.berkeley.edu/distmesh/)

% Reference: For derivation and fomula of the wave drag - 
%            Nycander, J. (2005). Generation of internal waves in the deep
%            ocean by tides. Journal of Geophysical Research, 
%            110(C10), C10028. http://doi.org/10.1029/2004JC002487

% Set some constants
% Earth Radius
R = 6378206.4;
% Coriolis coefficient
psi = 2*7.29212d-5;
% beta coefficient
beta = 1.455;

% Some arbitrary settings for the numerics 
rcutfac = 5 ; % Elements within |r - r'| < rcutfac*a will be used in the calculation of dJ/dx
acut = 2   ; % if a < acut (in km), ignore 

% Setting up parallel
if parallel > 0
    delete(gcp('nocreate'));
    pargroup = parpool(parallel) ; Pool = gcp('nocreate'); 
end
    
%% Computing the cutoff length of green's function
f = psi*sind(VX(:,2)) ; % Coriolis coefficients
dom = pi*sqrt(omega^2 - f.^2);
a = (beta./dom).*B.*Nm ;      % the cutoff length

%% Doing the projection
% Setting up projection
m_proj(proj,'lon',[ min(VX(:,1)) max(VX(:,1))],...
            'lat',[ min(VX(:,2)) max(VX(:,2))])
% Doing conversion
[XX(:,1),XX(:,2)] = m_ll2xy(VX(:,1),VX(:,2));      
XX = R * XX; % Convert to actual distances 
np = length(XX) ;  % Number of points

xc = (XX(EToV(:,1),:) + XX(EToV(:,2),:) + XX(EToV(:,3),:))/3 ; % centroid of each element

%% Setting up the quadrature and performing on the mesh
% Linear element
p = 1 ; % Do not change this value 
RefObj = InitRefOBJElVal( p ) ;

% 2-D quadrature
pcub = pcub_c*p ; % pcub_c order quadrature
RefObj2dCub = initvolcubintg( pcub, p, RefObj ) ;

% Find
% - Quadrature points for each elements
% - Bathymetric height (-B) at the quarature points
[hcub,xcub,ycub,jc] = getquadravalp1( RefObj2dCub, EToV, -B, XX ) ;

% Initialising the arrays
if length(hcut) > 1
    idv = find( (B > hcut(1)) & (B < hcut(2)) & ((a/1000) > acut) ) ;
else
    idv = find( (B > hcut)  & ((a/1000) > acut) ) ;
end
    
nj = length(idv) ;

dJ = zeros(np,2) ;
J  = zeros(np,1)  ;

%% Doing the computation (in parallel or serial)
if parallel > 0
    % Parallel
    idx = 0;
    num_p = 10*parallel;
    for jj = 1:10 %this for loop is for memory issues
        % Now lets evaluate the get_J function
        for ii = 1:parallel
            idx = idx + 1;
            ns = int64((idx-1)*nj/num_p)+1;
            ne = int64(idx*nj/num_p);
            fJ(idx) = parfeval(Pool,@get_J,2,XX(idv(ns:ne),:),-B(idv(ns:ne)),...
                a(idv(ns:ne)),xc,xcub,ycub,hcub,jc,rcutfac,RefObj2dCub);
        end
        disp('fetching the local evaluation of J and dJ')
        for ii = 1:parallel
            [idx_t, J_t, dJ_t] = fetchNext(fJ); % Get results into a cell array
            ns = int64((idx_t-1)*nj/num_p)+1;
            ne = int64(idx_t*nj/num_p);
            J(idv(ns:ne)) = J_t;
            dJ(idv(ns:ne),:) = dJ_t;
            disp(['fetched J and dJ for idx = ' num2str(idx_t)])
        end
        clear dJ
    end
    delete(pargroup);
else
    % Serial
    JvTemp   = zeros(nj,1) ; 
    dJvxTemp = zeros(nj,1) ;
    dJvyTemp = zeros(nj,1) ;
    tic 
    for i = 1: nj
        %
        iv = idv(i) ;
        xv = XX(iv,:) ;

        % 
        amv = a(iv) ;

        ielm = find(-dcircle( xc, xv(1), xv(2), rcutfac*amv ) > 0) ;

        if mod(i,ceil(nj/100)) == 0
            fprintf('node = %d ; nsum = %d ; amv = %f (km) \n', ...
                     iv, length(ielm), amv/1000 ) ;
        end

        % Bottom topography (increasing upward)
        bs = -B(iv) ;
        hs = hcub(:,ielm) - bs ; % I would use: hs = hcub(:,ielm)
        xs = xv(1) - xcub(:,ielm) ; % x - x'
        ys = xv(2) - ycub(:,ielm) ; % y - y' 

        % radius
        rrv = sqrt(xs.^2 + ys.^2) ;
        rrvinv = 1./rrv ; % 1/r

        dgar = diffgrnfunc( rrv, amv ) ; % dg_{a}/dr

        % dJ/dx
        vintg = RefObj2dCub.wc'*(hs.*dgar.*xs.*rrvinv) ;
        dJvxTemp(i) = vintg*jc(ielm) ;

        % dJ/dy
        vintg = RefObj2dCub.wc'*(hs.*dgar.*ys.*rrvinv) ;
        dJvyTemp(i) = vintg*jc(ielm) ;

        % J
        dgar = grnfunc( rrv, amv ) ;
        vintg = RefObj2dCub.wc'*(hs.*dgar) ;
        JvTemp(i) = vintg*jc(ielm) ;
    end
    % Enter into the global arrays
    dJ(idv,1) = dJvxTemp ;
    dJ(idv,2) = dJvyTemp ; 
    J(idv) = JvTemp ; 
    toc 
end
%EOF
end

function [J,dJ] = get_J(XX,h,a,xc,xcub,ycub,hcub,jc,rcutfac,RefObj2dCub)
    nj = length(h);
    J = zeros(nj,1);
    dJ = zeros(nj,2);
    for i = 1: nj
        %
        xv = XX(i,:);
        amv = a(i) ;

        ielm = find(-dcircle( xc, xv(1), xv(2), rcutfac*amv ) > 0) ;

        % Bottom topography (increasing upward)
        bs = h(i) ;
        hs = hcub(:,ielm) - bs ; % I would use: hs = hcub(:,ielm)
        xs = xv(1) - xcub(:,ielm) ; % x - x'
        ys = xv(2) - ycub(:,ielm) ; % y - y' 

        % radius
        rrv = sqrt(xs.^2 + ys.^2) ;
        rrvinv = 1./rrv ; % 1/r

        dgar = diffgrnfunc( rrv, amv ) ; % dg_{a}/dr

        % dJ/dx
        vintg = RefObj2dCub.wc'*(hs.*dgar.*xs.*rrvinv) ;
        dJ(i,1) = vintg*jc(ielm) ;
        

        % dJ/dy
        vintg = RefObj2dCub.wc'*(hs.*dgar.*ys.*rrvinv) ;
        dJ(i,2) = vintg*jc(ielm) ;

        % J
        dgar = grnfunc( rrv, amv ) ;
        vintg = RefObj2dCub.wc'*(hs.*dgar) ;
        J(i) = vintg*jc(ielm) ;
    end
end

