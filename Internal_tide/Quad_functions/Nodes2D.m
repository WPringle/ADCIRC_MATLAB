function [xx,yy,x,y] = Nodes2D( N )
% function [x,y,L1,L2,L3,warpf1,warpf2,warpf3,warp1,warp2,warp3] = Nodes2D( N )
% Taken from JS
%
alpcoef = [0.0000 0.0000 1.4152 0.1001 0.2751 0.9800 1.0999 ...
           1.2832 1.3648 1.4773 1.4959 1.5743 1.5770 1.6223 1.6258] ; 

if ( N < 16 )
    alpha = alpcoef(N) ;
else
    alpha = 5/3 ;
end
Np = (N + 1)*(N + 2)/2 ; 

L1 = zeros(Np,1) ;
L2 = zeros(Np,1) ; 
L3 = zeros(Np,1) ; 
sk = 1 ; 
for n = 1: N + 1
    for m = 1: N + 2 - n
        L1(sk) = (n - 1)/N ; 
        L3(sk) = (m - 1)/N ;
        sk = sk + 1 ; 
    end
end
% Generate uniform node distribution
L2 = 1.0 - L1 - L3 ; 
x = L3 - L2 ; 
y = (2*L1  -L3 - L2)/sqrt(3) ;

%
blend1 = 4*L2.*L3 ; 
blend2 = 4*L1.*L3 ;
blend3 = 4*L1.*L2 ;

warpf1 = Warpfactor(N, L3 - L2) ; 
warpf2 = Warpfactor(N, L1 - L3) ; 
warpf3 = Warpfactor(N, L2 - L1) ; 

warp1 = blend1.*warpf1.*(1 + (alpha*L1).^2) ;
warp2 = blend2.*warpf2.*(1 + (alpha*L2).^2) ;
warp3 = blend3.*warpf3.*(1 + (alpha*L3).^2) ;

xx = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3 ;
yy = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3 ; 