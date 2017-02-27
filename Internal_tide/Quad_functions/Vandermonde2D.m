function [V2D] = Vandermonde2D( N, r, s)
%
% Given (r,s) find Vandermode matrix 
V2D = zeros(length(r),(N+1)*(N+2)/2) ;

sk = 1 ;
for i = 0: N 
    for j = 0: N - i
        for l = 1: length(r)
            V2D(l,sk) = Simplex2DP(r(l),s(l),i,j) ;
        end
        sk = sk + 1 ; 
    end
end