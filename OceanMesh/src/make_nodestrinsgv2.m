clearvars; close all;
% Make node strings for ocean and land boundaries only. 
% improved by kjr, UND, CHL, Sept 2017
outfilename = 'COMBINED3.14'; 
% Nodes for breaking outer polygon into ocean and mainland in ccw order
ob_se =[3450223  2336406 ];
ob_se =[];

disp(['Reading in mesh ',outfilename,'...'])
tic
[ev,pv,B] = readfort14(outfilename,0);
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

%  For ocean & mainland boundary
nope = 0; neta = 0;
if(isempty(ob_se))
    % If no specified ocean start and end points, then everything is mainland.
    nbou = nbou + 1;
    nvell(nbou) = length(vso);
    nvel = nvel + nvell(nbou);
    nbvv(1:nvell(nbou),nbou) = idv';
    ibtype(nbou) = 20;
    plot(vso(:,1),vso(:,2),'k-'); 
else
    for ob = ob_se'
        % Find start and end of ocean boundaries
        idx1 = find(idv == ob(1));
        idx2= find(idv == ob(2));
        
        if(idx1 < idx2)
            I_s=idx1;
            I_e=idx2;
        else
            I_s=idx2;
            I_e=idx1;
        end
        
        if ~isempty(I_s) && ~isempty(I_e)
            nope = nope + 1;
            nvdll(nope) = length(idv(I_s:I_e));
            neta = neta + nvdll(nope);
            ibtypee(nope) = 0;
            nbdv(1:nvdll(nope),nope) = idv(I_s:I_e)';
            plot(pv(nbdv(1:nvdll(nope),nope),1),pv(nbdv(1:nvdll(nope),nope),2))
            
            % Get mainland boundary
            nbou = nbou + 1;
            ibtype(nbou) = 20;
            I_s = I_e;
            if nope < size(ob_se,1)
                I_e = find(idv == ob_se(nope+1,1));
            else
                I_e = find(idv == ob_se(1,1));
            end
            if I_e > I_s
                nvell(nbou) = length(idv(I_s:I_e));
                nbvv(1:nvell(nbou),nbou) = idv(I_s:I_e)';
            else
                nvell(nbou) = length(idv(I_s:end)) + length(idv(1:I_e));
                nbvv(1:nvell(nbou),nbou) = [idv(I_s:end)'; idv(1:I_e)'];
            end
            nvel = nvel + nvell(nbou);
            hold on ;plot(pv(nbvv(1:nvell(nbou),nbou),1),pv(nbvv(1:nvell(nbou),nbou),2))
        end
        
    end   
 end
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
writefort14( [outfilename,'wNS.grd'] , ev, pv, B, opedat , boudat ,'grid' ) ;