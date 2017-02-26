function [fcub,xcub,ycub,jx] = getquadravalp1( RefObj2dCub, EToV, fxy, XX )
%
% For linear element only
%
% Given  fxy at FE nodes return fxy at quadrature nodes 
%    dim(fxy) =  (np,1) 
%    dim(EToV) = (ne,3)
%
%    dim(fcub) = (npcub,ne)
% if the FE nodes, i.e. XX is given, the function also returns the
% quadrature nodes, and Jacobian 
%   dim(xcub) = (npcub,ne)
%   dim(ycub) = (npcub,ne)
%   dim(jx) = (ne,1)
%
%  DW:
%

if ( RefObj2dCub.p ~= 1 )
    disp('Error: RefObj2dCub.p ~= 1. This function is for linear elements only') ;
    fcub = [] ;
    xcub = [] ; 
    ycub = [] ; 
    return ; 
end

% 
ne = length(EToV) ; 
fcub = RefObj2dCub.Vcub*fxy(EToV') ;

xcub = [] ; 
ycub = [] ;
jx = [] ; 
if ( nargin == 4 )
  xcub = RefObj2dCub.Vcub*reshape(XX(EToV',1),3,ne) ;  
  ycub = RefObj2dCub.Vcub*reshape(XX(EToV',2),3,ne) ; 
  
  % Find Jacobian
  xr = 0.5*(XX(EToV(:,2),:) - XX(EToV(:,1),:)) ; % x_{r}, y_{r}
  xs = 0.5*(XX(EToV(:,3),:) - XX(EToV(:,1),:)) ; % x_{s}, y_{s}
  
  % Jacobian
  jx = xr(:,1).*xs(:,2) - xs(:,1).*xr(:,2) ; 
end