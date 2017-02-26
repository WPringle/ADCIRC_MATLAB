% find band width of the grid
f14 = 'prviv19h.grd.14' ; 
[etov,vx,~] = readfort14( f14 ) ;

%
[vtov, nnv,] = NodeConnect2D( etov ) ;

[ir,ic] = size(vtov) ;
icolv = ones(ir,1)*[1:ic] ;

idr = find(vtov > 0) ;
idrow = vtov(idr) ; 
idcol = icolv(idr) ; 
vtovmat = sparse( idrow, idcol, idrow*0 + 1) ; 

[lb,ub] = bandwidth( vtovmat ) ; 