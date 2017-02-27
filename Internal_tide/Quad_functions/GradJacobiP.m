function [dP] = GradJacobiP( r, alpha, beta, N)
% Taken from JSH

dP = zeros(length(r), 1) ;
if ( N == 0 )
    dP(:,:) = 0.0 ; 
else
    dP = sqrt(N*(N + alpha + beta + 1))*JacobiP(r(:),alpha+1,beta+1,N-1) ;
end

return ;