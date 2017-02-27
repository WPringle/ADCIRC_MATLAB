function RefObjCub = initvolcubintg( pcub, p, RefObj )
%
% Quadrature for volume integral evaluttion
%   pcub - order of quadrature 
%   p    - degree of interpolant
%   RefObj - Reference element matrices associated with interpolant degree
%
%  
%  xc -- Quadrature points
%  wc -- Quadrature weights
%  ncup -- number of quadrature points
%  Vcup -- Interpolation matrix, i.e. Vcup*f --> f_{c}
%  V2Dr -- Map function valies at the nodes to its derivative with
%          respect to r at the quadrature point  (nodal to nodal)  
%  V2Ds -- Map function valies at the nodes to its derivative with
%          respect to s at the quadrature points (nodal to nodal)   
%

cub2D = GetQuadrature2D( pcub ) ;
cubw = cub2D(:,3) ; 

% Vcub == V_{c} V^{-1}
Vcub = Vandermonde2Dc( cub2D, p, RefObj ) ;

% V_{,r}, V_{,s}
[V2Dr,V2Ds] = GradVandermonde2D( p, cub2D(:,1), cub2D(:,2) ) ;

np = find( RefObj.p == p ) ;

% V_{,r} = V_{,r} V^{-1} 
% V_{,s} = V_{,s} V^{-1}
V2Dr = V2Dr*RefObj.REl(np).InvRefVanderMat ;
V2Ds = V2Ds*RefObj.REl(np).InvRefVanderMat ;

%
% V_{,rr}, V_{,ss}, V_{,rs}
%
%  (V_{**})_{ij} = \phi_{j,**}(x_{c,i})
[V2Drr,V2Dss,V2Drs] = DGradVandermonde2D( p, cub2D(:,1), cub2D(:,2) ) ;
V2Drr = V2Drr*RefObj.REl(np).InvRefVanderMat ;
V2Dss = V2Dss*RefObj.REl(np).InvRefVanderMat ;
V2Drs = V2Drs*RefObj.REl(np).InvRefVanderMat ;


% Quadratures
RefObjCub.pcub = pcub ; 
RefObjCub.xc = cub2D(:,1:2) ; 
RefObjCub.wc = cubw ; 
RefObjCub.ncub = length(cubw) ;

% For calcuation of integral term
% using a quadrature rule
RefObjCub.p = p ;
RefObjCub.Vcub = Vcub ;
RefObjCub.VDc(1).Vd = V2Dr ;
RefObjCub.VDc(2).Vd = V2Ds ;

%
% For a calculation involving with the second derivative 
% of the basis functions
%
% ---- Added on Dec 1, 2014
%
RefObjCub.V2Dc(1,1).Vd = V2Drr ;
RefObjCub.V2Dc(1,2).Vd = V2Drs ;
RefObjCub.V2Dc(2,1).Vd = V2Drs ; 
RefObjCub.V2Dc(2,2).Vd = V2Dss ;
