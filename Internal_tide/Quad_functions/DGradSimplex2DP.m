function [dmodedrr,dmodedss,dmodedrs] = DGradSimplex2DP(r, s, id, jd)
% 
%  Second derivative:
%
%    \phi_{r,r}(r,s), \phi_{r,s}(r,s), \phi_{s,s}(r,s)
% 
%  DW: 
%  Mental note
%
%     - Singularity must be treated carefully
%
%  
n = length(r) ;
a = zeros(n,1) ;

eps = 1e-14 ; 
for ii = 1: n
    if ( abs(s(ii) - 1) > eps )
        a(ii) = 2*(1 + r(ii))/(1 - s(ii)) -  1 ;
    else
        a(ii) = -1 ;
    end   
end
b = s ;

fa = JacobiP( a, 0, 0, id ) ; dfa = GradJacobiP( a, 0, 0, id ) ; 
gb = JacobiP( b, 2*id+1, 0, jd ) ; dgb = GradJacobiP( b, 2*id+1, 0, jd ) ; 

% 
d2fa = DnGradJacobiP( a, 2, 0, 0, id ) ; 
d2gb = DnGradJacobiP( b, 2, 2*id + 1, 0, jd ) ; 

%
% phi_{r,r}
dmodedrr = gb.*d2fa ; 
if ( id > 1 ) 
   dmodedrr = 2^(5/2)*(dmodedrr.*(1 - b).^(id - 2)) ;
end
% 

%
% phi_{r,s}
%
dmodedrs = 0*gb ;
if ( id > 1 )
    dmodedrs = gb.*d2fa.*(0.5*(1 + a)) ;
    
    dmodedrs = 2^(5/2)*(dmodedrs.*(1 - b).^(id - 2)) ;
end

tmp = dfa.*dgb ; 
if ( id > 0 )
    %
    tmp = 2^(3/2)*(tmp.*(1 - b).^(id - 1)) ;
end
dmodedrs = dmodedrs + tmp ;

if ( id > 1 )
    tmp = dfa.*gb ;
    
    tmp = -2^(3/2)*((id - 1)*tmp.*(1 - b).^(id -2)) ;
    
    dmodedrs = dmodedrs + tmp ;
end

%
%
% phi_{ss}
%
dmodedss = gb*0 ;
if ( id > 1 )
    dmodedss = d2fa.*gb.*(0.25*(1 + a).^2) ;
    %
    
    dmodedss = 4*dmodedss.*(1 - b).^(id - 2) ;
end

% + 
if ( id > 0 ) 
   tmp = dfa.*dgb.*(0.5*(1 + a)) ; 

   tmp = (2 + 2)*((1 - b).^(id - 1)).*tmp ;
   
   dmodedss = dmodedss + tmp ;
end

%
if ( id > 1 )
    tmp = dfa.*gb.*(0.5*(1 + a)) ;
    tmp = (-2*(id - 2) - 2*id)*tmp.*(1 - b).^(id - 2) ;
    
    dmodedss = dmodedss + tmp ;
end

% 
if ( id > 0 )
    %
    % tmp = fa.*dgb ;
    %
    tmp = (-id - id)*fa.*dgb.*(1 - b).^(id - 1) ;
    
    dmodedss = dmodedss + tmp ;
end

if ( id > 1 )
    %
    tmp = fa.*gb ;
    tmp = id*(id - 1)*tmp.*(1 - b).^(id - 2) ;
    
    dmodedss = dmodedss + tmp ;
end

tmp = fa.*d2gb ;
if ( id > 0 ) 
    tmp = tmp.*(1 - b).^(id) ;
end
dmodedss = dmodedss + tmp ;

% D_{ss}
dmodedss = sqrt(2)*dmodedss ; 

return ;