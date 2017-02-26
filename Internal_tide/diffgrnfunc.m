function  gar = diffgrnfunc( r, a )
%
% Derivative of scaled Green function
%
%  dg_{a}/dr; g_{a} = (1/a)*G(r/a)
%
%  G = 1/x - (\pi/2)*exp(-x^2/8)I_{0}(x^{2}/8)
%
% From Maple:
%  dG/dx = -0.1e1 / x ^ 2 + sqrt(pi) * x * exp(-x ^ 2 / 0.8e1) * besseli(0.0e0, x ^ 2 / 0.8e1) / 0.8e1 ...
%           - sqrt(pi) * exp(-x ^ 2 / 0.8e1) * besseli(0.1e1, x ^ 2 / 0.8e1) * x / 0.8e1
%
aiv = (1/a) ;
x  = aiv*r ;
x2 = (x.*x)/8 ;

fcv = sqrt(pi)/8 ;
gar = -1./(x.*x) + fcv*x.*exp(-x2).*besseli(0,x2) - fcv*x.*exp(-x2).*besseli(1,x2) ;

%
% d g_{a}/dr 
%
gar = (aiv*aiv)*gar ;

end