function [dmodedr,dmodeds] = GradSimplex2DP(r, s, id, jd)
%
% Taken from JSH, differentation of modal basis
%       d \psi_{(id,jd)}(r,s)/dr     0 <= id <= N, (id+jd) <= N 
%
%
n = length(r) ;
a = zeros(n,1) ; 

for ii = 1: n
    if ( abs(s(ii) - 1) > 1e-14 )
        a(ii) = 2*(1 + r(ii))/(1 - s(ii)) -  1 ;
    else
        a(ii) = -1 ;
    end
end
b = s ; 

fa = JacobiP( a, 0, 0, id ) ; dfa = GradJacobiP( a, 0, 0, id ) ; 

gb = JacobiP( b, 2*id+1, 0, jd ) ; dgb = GradJacobiP( b, 2*id+1, 0, jd ) ; 

dmodedr = dfa.*gb ; 
if ( id > 0 )
    dmodedr = dmodedr.*(2*(1 - b).^(id - 1)) ;
    % dmodedr = dmodedr.*((0.5*(1 - b)).^(id - 1)) ;
end

dmodeds = dfa.*(gb.*(1 + a)) ;
% dmodeds = dfa.*(gb.*(0.5*(1 + a))) ;
if ( id > 0 )
    dmodeds = dmodeds.*(1-b).^(id - 1) ;
    % dmodeds = dmodeds.*((0.5*(1-b)).^(id - 1)) ;
end

tmp = dgb.*((1 - b).^id) ;
% tmp = dgb.*((0.5*(1 - b)).^id) ;
if ( id > 0 )
     tmp = tmp - id*gb.*((1 - b).^(id - 1)) ;
     % tmp = tmp - 0.5*id*gb.*((0.5*(1 - b)).^(id - 1)) ;
end
dmodeds = dmodeds + fa.*tmp ;

dmodedr = sqrt(2)*dmodedr ;
dmodeds = sqrt(2)*dmodeds ;

return ;