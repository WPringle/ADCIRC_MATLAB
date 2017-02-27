function [x] = JacobiGL(alpha, beta, N)
%
x = zeros(N+1,1) ;
if ( N == 1 )
    x(1) = -1.0 ; 
    x(2) = 1.0 ;
    return ;
end

[xint,w] = JacobiGQ(alpha+1,beta+1,N-2) ;
x = [-1, xint', 1]' ; 

return ;
