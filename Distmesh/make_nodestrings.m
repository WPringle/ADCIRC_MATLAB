% Load mesh
clearvars; close all;

[ev,pv,B,~,~,~] = readfort14('INDPAC_sponge.grd');

% Nodes for breaking outer polygon into ocean and mainland
ob_se = [3256547 3571732;
        258906 1460086;
        2894522 2838695];

% Extrace edges on the boundaries from a given mesh, pv, ev %
[etbv,vxe,etoe,etof] = extdom_edges( ev, pv ) ;
%subplot(2,2,1) ; triplot( ev, pv(:,1), pv(:,2) ) ;

if ~isempty(ob_se)
    iedbeg  = knnsearch(pv(etbv(1,:),:),pv(ob_se(1,1),:));
else
    iedbeg  = 1 ;
end
    
start = 1;
nbou = 0;
nvel = 0;
% Extract the outter domain 
while ~isempty(etbv)
    pacw = 0; nn = 0;
    if start == 0
        iedbeg = 1;
    end
    while pacw == 0
        nn = nn + 1;
        if nn == 1
            ipsbeg = 1 ; % polygon travese in 'ipsbeg' --> 'ipsend' 
            ipsend = 2 ; %
        else
            ipsbeg = 2 ; % polygon travese in 'ipsbeg' --> 'ipsend' 
            ipsend = 1 ; %   
        end
        [vso,idv,ide] = extdom_polygon( etbv, pv, iedbeg, ipsbeg, ipsend ) ;
        if length(vso) < 6
        end
        pacw = ~ispolycw(vso(:,1),vso(:,2));
        if start == 1 && pacw == 0
            break; 
        end
        % vso(:,2) - coordinate of an extracted polygon 
        % idv  -- indices of extracted ploygon in the global mesh,  vso = pv(idv,:)  
        % ide  -- indeces of edges in etbv that constitute vso,  
        % ie. vso = vxe(reshape(etbv(:,ide),length(ide)*2),:)  
    end
    plot(vso(:,1),vso(:,2))
    hold on
    if start == 1
        if ~isempty(ob_se)
            % ocean & mainland
            nope = 0; neta = 0;
            for ob = ob_se'
                % Get ocean
                I_s = find(idv == ob(1));
                I_e = find(idv == ob(2));
                nope = nope + 1;
                nvdll(nope) = length(idv(I_s:I_e));
                neta = neta + nvdll(nope);
                ibtypee(nope) = 0;
                nbdv(1:nvdll(nope),nope) = idv(I_s:I_e)'; 
                plot(pv(nbdv(1:nvdll(nope),nope),1),pv(nbdv(1:nvdll(nope),nope),2))
                
                % Get mainland
                nbou = nbou + 1;
                ibtype(nbou) = 20;
                I_s = I_e;
                if nope < length(ob_se)
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
                plot(pv(nbvv(1:nvell(nbou),nbou),1),pv(nbvv(1:nvell(nbou),nbou),2))
            end
        else
            % mainland
            nbou = nbou + 1;
            nvell(nbou) = length(vso);
            nvel = nvel + nvell(nbou);
            nbvv(1:nvell(nbou),nbou) = idv'; 
            ibtype(nbou) = 20;
        end
        start = 0;
    else
        % island
        nbou = nbou + 1;
        nvell(nbou) = length(vso);
        nvel = nvel + nvell(nbou);
        nbvv(1:nvell(nbou),nbou) = idv'; 
        ibtype(nbou) = 21;
    end
    etbv(:,ide) = []; 
end

if ~isempty(ob_se)
    % ocean boundary
    opedat.nope = nope ; 
    opedat.neta = neta ;
    opedat.nvdll = nvdll ;
    opedat.ibtypee = ibtypee ;
    opedat.nbdv = nbdv ; 
else
    opedat = [];
end


% land boundary
boudat.nbou = nbou ;
boudat.nvel = nvel ;
boudat.nvell = nvell ;
boudat.ibtype = ibtype ;
boudat.nbvv = nbvv ;  

writefort14( 'INDPAC_sponge_ns.grd' , ev, pv, ...
              B, opedat , boudat ,'grid' ) ;
