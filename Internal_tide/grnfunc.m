function  gar = grnfunc( r, a )
%
% Scaled Green function
%
%  g_{a} = (1/a)*G(r/a)
%
%  G = 1/x - (\pi/2)*exp(-x^2/8)I_{0}(x^{2}/8)
%

aiv = (1/a) ; 
x = aiv*r ;
x2 = (x.*x)/8 ; % (1/8)*(r/a)^2

fcv = 0.5*sqrt(pi) ;
gar = -fcv*exp(-x2).*besseli(0,x2) ; % bsesseli -- modify Bessel of the first kind

gar = aiv*(gar + (1./x)) ;


