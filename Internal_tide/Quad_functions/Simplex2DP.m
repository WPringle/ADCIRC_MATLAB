function [P] = Simplex2DP(r, s, i, j)
%
% Taken from JS Hesthaven
%
%  Calculate value of 
%
%    \sqrt(2) P^{i} P^{2i + 1,j} (1 - b)^{i}
%
n = length(r) ;
a = zeros(n,1) ;

for  ii = 1: n
    if ( abs(s(ii) - 1) > 1e-14 )
        a(ii) = 2*(1 + r(ii))/(1 - s(ii)) -  1 ;
    else
        a(ii) = -1 ;
    end
end

b = s ; 

Pi = JacobiP(a, 0, 0, i) ;
Pj = JacobiP(b, 2*i + 1, 0, j) ;

P = sqrt(2)*Pi.*Pj.*(1 - b).^i ;



