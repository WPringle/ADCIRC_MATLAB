function [V1D] = Vandermonde1D( N, r ) 
% Taken from JSH 

V1D = zeros(length(r),N+1) ;
for j = 1:N +1 
    V1D(:,j) = JacobiP(r(:), 0, 0, j -1) ;
end

return ;  
