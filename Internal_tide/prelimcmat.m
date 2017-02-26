clear ;

load('Nb_Nm_values.mat') ;
whos ;

% Earth Radius
R = 6378206.4;

% Coriolis coefficient
psi = 2*7.29212d-5;

% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
tidefreq = 12.4206012 ; % Hr
omega = 2*pi/(tidefreq*3600) ; % M2 tidal frequency

bet = 1.455 ;
f = psi*sin(pi*VX(:,2)/180) ; % Coriolis coefficients

dom = pi*sqrt(omega^2 - f.^2) ;
a = (bet./dom).*B.*Nm ; % radius

%
% -- CPP projection (May not be the best idea ??)
% -- (lon0,lat0) - in rad 
XX = VX*pi/180 ;
XX(:,1) = R * (XX(:,1) - lon0) * cos(lat0);
XX(:,2) = R * XX(:,2) ;
np = length(XX) ;  % Number of points

%
% [amax,ib] = max( a ) ;
% rr = sqrt((XX(:,1) - XX(ib,1)).^2 + (XX(:,2) - XX(ib,2)).^2) ;
%
% gar = grnfunc( rr, amax ) ;    
% dgar = diffgrnfunc( rr, amax ) ;
%
%

% Linear element
p = 1 ; % Do not change this value 
RefObj = InitRefOBJElVal( p ) ;

% 2-D quadrature
pcub = 2*p ; % 3 - point quadrature
RefObj2dCub = initvolcubintg( pcub, p, RefObj ) ;

% Find
% - Quadrature points for each elements
% - Bathymetric depth at the quarature points
[bcub,xcub,ycub,jc] = getquadravalp1( RefObj2dCub, EToV, B, XX ) ;

rcutfac = 2.5 ; % Elements within |r - r'| < rcutfac*a will be used in the calculation of dJ/dx
hcut = 140 ;     % Find dJ/dx, dJ/dy, J when H > hcut
acut = 2   ; % if a < acut (in km), ignore 
xc = (XX(EToV(:,1),:) + XX(EToV(:,2),:) + XX(EToV(:,3),:))/3 ; % centroid of each element

idv = find( (B > hcut)  & ((a/1000) > acut) ) ;

nj = length(idv) ;


dJv = zeros(2,np) ;
Jv = zeros(1,np) ;

JvTemp = zeros(1,nj) ; 
dJvxTemp = zeros(1,nj) ;
dJvyTemp = zeros(1,nj) ;


% pargroup = parpool(4) ;
tic 

% parfor i = 1: nj
for i = 1: nj
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
dJv(1,idv) = dJvxTemp ;
dJv(2,idv) = dJvyTemp ; 
Jv(idv) = JvTemp ; 
toc 

%
% delete(pargroup);  


