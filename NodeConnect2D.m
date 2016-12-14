function [ vtov, nnv, vtoe, nn, etab ] = NodeConnect2D( etov )
%
% Table listing element sharing the node
%
nel = length(etov) ; 
nv = max(max(etov)) ;

icc = (1:nel)'*ones(1,3) ;
etab = sparse(icc, etov, etov*0 + 1 ) ;

nn = full(sum(etab)) ;

tic ;
nme = max(nn) ; 
vtoe = zeros(nme,nv) ;
for i = 1: nv
   vtoe(1:nn(i),i) = find(etab(:,i)) ;    
end
toc ;

tic ;
% node connectivity
vtov = zeros(nme+3,nv) ;

ev = etov' ;

nnv = zeros(nv,1) ; 
for i = 1: nv
   ineigh = vtoe(1:nn(i),i) ;
    
   val = unique(ev(:,ineigh)) ;
   
   nnv(i) = length(val) ; 
   
   vtov(1:nnv(i),i) = val ;
end
toc ;

end

