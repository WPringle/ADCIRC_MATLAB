function  Vc = Vandermonde2Dc( CubQ2D, p, RefEl ) 
%
% Get Van Der Monde matrix for associate with 2D Quadratures 
%
% Output:
%  Vc - interpolate function at the quadrature points from
%       nodal values associate with a degree-p interpolant, i.e.  
%       Vc*f ---> fc
%
[N,ic] = size( CubQ2D ) ; 
r = CubQ2D(:,1) ; 
s = CubQ2D(:,2) ; 

V2D = zeros(N,(p+1)*(p+2)/2) ;
Vc  = zeros(N,(p+1)*(p+2)/2) ; 

sk = 1 ;
for i = 0: p 
    for j = 0: p - i
        %
        for l = 1: N
            V2D(l,sk) = Simplex2DP(r(l),s(l),i,j) ;
        end
        sk = sk + 1 ;
    end
end

% Get interpolation matrix %
np = find( RefEl.p == p ) ; 
Vc = V2D*RefEl.REl(np).InvRefVanderMat ; 