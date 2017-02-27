function [J,dJ] = Compute_J_Nycander(EToV,VX,B,Nm,omega,...
                                     beta,pcub_c,hcut,proj,parallel)
% Computes J, and its derivatives for use in Nycander (2005) wave drag
% formula on the finite-element grid using quadrature of arbitrary degree
% function [J,dJ] = Compute_J_Nycander(EToV,VX,B,Nb,Nm,beta,omega,hcut,proj,parallel)   
% 
% Inputs:  EToV - the element table (Num_ele x 3)     
%          VX   - the vertex positions (Num_node x 2)
%          B    - the bathymetric depth at each node (Num_node x 1)
%          Nm   - the depth-averaged Brunt–Väisälä buoyancy frequency at
%                 each node (Num_node x 1)
%         omega - the circular frequency (rad/s) of the tidal component of
%                 interest (usually the M2 one)
%          beta - a numerical coefficient to determine the cutoff length of
%                 the Green's function kernel (beta = 1.455 is typically
%                 used as defined by Nycander, 2005)
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
%            - Distmesh toolbox for dcircle function

% Reference: For derivation and fomula of the wave drag - 
%            Nycander, J. (2005). Generation of internal waves in the deep
%            ocean by tides. Journal of Geophysical Research, 
%            110(C10), C10028. http://doi.org/10.1029/2004JC002487

% Set some constants
% Earth Radius
R = 6378206.4;
% Coriolis coefficient
psi = 2*7.29212d-5;

% Some arbitrary settings for the numerics 
rcutfac = 2.5 ; % Elements within |r - r'| < rcutfac*a will be used in the calculation of dJ/dx
acut = 2   ; % if a < acut (in km), ignore 

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
pcub = pcub_c*p ; % 3 - point quadrature
RefObj2dCub = initvolcubintg( pcub, p, RefObj ) ;

% Find
% - Quadrature points for each elements
% - Bathymetric depth at the quarature points
[bcub,xcub,ycub,jc] = getquadravalp1( RefObj2dCub, EToV, B, XX ) ;

%% Initialising the arrays and opening parallel pool if necessary
idv = find( (B > hcut)  & ((a/1000) > acut) ) ;

nj = length(idv) ;

dJ = zeros(np,2) ;
J  = zeros(np,1)  ;

JvTemp   = zeros(nj,1) ; 
dJvxTemp = zeros(nj,1) ;
dJvyTemp = zeros(nj,1) ;

if parallel > 0
    pargroup = parpool(parallel) ;
end

%% Perform the computation looping over all necessary nodes
tic 
parfor i = 1: nj
%for i = 1: nj
    %
    iv = idv(i) ;
    xv = XX(iv,:) ;
   
    % 
    amv = a(iv) ;
   
    ielm = find(-dcircle( xc, xv(1), xv(2), rcutfac*amv ) > 0) ;

    fprintf('node = %d ; nsum = %d ; amv = %f (km) \n', iv, length(ielm), amv/1000 ) ;

    % Bottom topography (increasing upward)
    bs = -B(iv) ;
    hs = -bcub(:,ielm) - bs ; % I would use: hs = -bcub(:,ielm)
    xs =  xv(1) - xcub(:,ielm) ; % x - x'
    ys =  xv(2) - ycub(:,ielm) ; % y - y' 

    % radius
    rrv = sqrt(xs.^2 + ys.^2) ;
    rrvinv = 1./rrv ; % 1/r

    dgar = diffgrnfunc( rrv, amv ) ; % dg_{a}/dr

    % dJ/dx
    vintg = RefObj2dCub.wc'*(hs.*dgar.*xs.*rrvinv) ;
    % dJv(1,iv) = vintg*jc(ielm) ;
    dJvxTemp(i) = vintg*jc(ielm) ;
    
    % dJ/dy
    vintg = RefObj2dCub.wc'*(hs.*dgar.*ys.*rrvinv) ;
    % dJv(2,iv) = vintg*jc(ielm) ;
    dJvyTemp(i) = vintg*jc(ielm) ;

    % J
    dgar = grnfunc( rrv, amv ) ;
    vintg = RefObj2dCub.wc'*(hs.*dgar) ;
    %Jv(iv) = vintg*jc(ielm) ;
    JvTemp(i) = vintg*jc(ielm) ;
end

%% Enter into the global arrays
dJ(idv,1) = dJvxTemp ;
dJ(idv,2) = dJvyTemp ; 
J(idv) = JvTemp ; 
toc 
if parallel > 0 
    delete(pargroup);
end
%EOF
end


