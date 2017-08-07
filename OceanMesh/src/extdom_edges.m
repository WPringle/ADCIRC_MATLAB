function  [etbv,vxe,etoe,etof] = extdom_edges( ev, pv )
%
% Given Element table, extract edde 
%

[etoe, etof] = Connect2D( ev ) ;

% Extract the boundires
geid = find(~abs(etoe - [1:length(etoe)]'*ones(1,3))') ;

edelm = ceil(geid/3) ;

edid = mod(geid,3)   ;
edid(edid == 0) = 3  ;

edge = [ 1 2 3
         2 3 1 ] ;

ned = length(edid) ;

%
% Gather edge table
etabl = edge(:,edid)  ;
etbv  = zeros(size(etabl)) ;
for i = 1: ned
    etbv(:,i) = ev(  edelm(i), etabl(:,i) ) ;
end

idx = unique(etbv(:)) ;
vxe = pv(idx,:) ;

end
