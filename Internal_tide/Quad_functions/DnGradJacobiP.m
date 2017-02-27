function [dP] = DnGradJacobiP( r, dn, alpha, beta, N)
% 
%  d^(n) J^{alpha,beta}_{N}(r) /dr^{n}
% 
% DW
%
dP = zeros(length(r), 1) ;
if ( N <= 0  )
    dP(:) = 0.0 ;
    return ;
end

if ( dn == 1 ) 
    dP = GradJacobiP( r(:), alpha, beta, N ) ;
    
    return ;
end

dP = sqrt( N*(N + alpha + beta + 1) )*...
    DnGradJacobiP( r, dn - 1, alpha + 1, beta + 1, N - 1 ) ;   

% else
%    % dP(:) = 1.0 ;
%    tmp = 1.0 ;
%    
%    %for dn = 1: dn - 1
%    %    ii = dn - 1 ;
%    %    %
%    %    tmp = sqrt((N - ii)*((N - ii) + (alpha + ii) + (beta + ii) + 1))*tmp ; 
%    % 
%    %end
%    %dP = GradJacobiP( r(:), alpha + dn - 1, beta + dn - 1, N - dn ) ; 
%    
% end

return ;