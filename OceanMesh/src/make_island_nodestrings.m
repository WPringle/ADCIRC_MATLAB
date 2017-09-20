clearvars; close all;
% Make node strings for ocean and land boundaries only. 
% improved by kjr, UND, CHL, Sept 2017
outfilename = 'fixed.grd';
% Nodes for breaking outer polygon into ocean and mainland in ccw order
ob_se =[];

disp(['Reading in mesh ',outfilename,'...'])
tic
[ev,pv,B] = readfort14([outfilename]);
toc

disp('Tracing boundary in winding order...')
[etbv,vxe]              = extdom_edges2(ev,pv); 
tic
[poly,poly_idx,max_ind] = extdom_polygon(etbv,pv,1);
toc

disp('Identifying and storing boundaries...')
nbou = 0;
nvel = 0;
start= 1;
% the largest polygon will be a combination of ocean and mainland
% boundaries. Deal with this first, then remove it from the polygon
idv = poly_idx{max_ind};
vso = poly{max_ind};

poly(max_ind) = [];
poly_idx(max_ind)=[];

% loop through the remaining polygons
for poly_count = 1 : length(poly)
    vso = poly{poly_count};
    idv = poly_idx{poly_count};
    % islands
    nbou = nbou + 1;
    nvell(nbou) = length(vso);
    nvel = nvel + nvell(nbou);
    nbvv(1:nvell(nbou),nbou) = idv';
    ibtype(nbou) = 21;
    hold on ;plot(pv(nbvv(1:nvell(nbou),nbou),1),pv(nbvv(1:nvell(nbou),nbou),2))
end

disp('Forming opendat and boudat...'); 
if ~isempty(ob_se)
    % ocean boundary
    opedat.nope = nope ;
    opedat.neta = neta ;
    opedat.nvdll = nvdll ;
    opedat.ibtypee = ibtypee ;
    opedat.nbdv = nbdv ;
else
    % ocean boundary
    opedat.nope = 0 ;
    opedat.neta = 0 ;
    opedat.nvdll = 0 ;
    opedat.ibtypee = 0 ;
    opedat.nbdv = 0;    
end


% land boundary
boudat.nbou = nbou ;
boudat.nvel = nvel ;
boudat.nvell = nvell ;
boudat.ibtype = ibtype ;
boudat.nbvv = nbvv ;  

disp('Writing fort.14 with nodestrings...'); 
writefort14( [outfilename,'wIslandNS.grd'] , ev, pv, B, opedat , boudat ,'grid' ) ;