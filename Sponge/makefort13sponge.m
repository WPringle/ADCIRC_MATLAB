function [sponge,opedat,boudat,pv,B] = makefort13sponge(f13file,f14file,...
                                      period,frac,spngtype,F,rat,opv,write)

g = 9.81;  % gravity
alpha = 2; % second order polynomial

%% Get the nodes for each sponge zone
[sponge,opedat,boudat,pv,B] = get_latlon_for_sponge_zone(f14file,...
                                                        period,frac,write);

%% make the sigma based on spongetype and coefficients
sigma = []; idspg_node = [];
for op = opv
    or_abs = abs(sponge(op).orientation);
    or_sign = sign(sponge(op).orientation);
    % the x (or y) distance from the open boundary
    [~,d] = knnsearch(...
     max(pv(opedat.nbdv(1:opedat.nvdll(op),op),or_abs))*(1-or_sign)/2 + ...
     min(pv(opedat.nbdv(1:opedat.nvdll(op),op),or_abs))*(1+or_sign)/2,...
                      sponge(op).pv(:,or_abs));
    % change orientation so origin is at the start of sponge
    d = max(sponge(op).L-d,0);
    % need to get metre equivalents of L & d (roughly)
    L = sponge(op).L*1d5; d = d*1d5;
    % spongetype polynomial or hyperbola
    if strcmp(spngtype,'poly')
        sigma_m = -sqrt(g*sponge(op).H)*(alpha+1)*log(1/F)/(L*rat^(alpha+1));  
        sigma_n = sigma_m*(d/L).^alpha;
    elseif  strcmp(spngtype,'hyper')
        sigma_m = log(1/F)/(-log(1-rat)-rat);
        sigma_n = sigma_m*sqrt(g*sponge(op).H)*d./(L*(L-d));
    end
    idspg_node = [idspg_node; sponge(op).idx];
    sigma = [sigma; sigma_n];
end

%% Set fort13 structure and write out to file
% User-defined input
f13dat.AGRID = 'spatial_attributes_description';
f13dat.NumOfNodes = length(pv);
f13dat.nAttr = 1; natb = f13dat.nAttr;
f13dat.userval.Atr(natb).AttrName = 'sponge_generator_layer' ; 
f13dat.userval.Atr(natb).usernumnodes = length(idspg_node) ;
f13dat.userval.Atr(natb).Val = [ idspg_node'; sigma'; (sigma*0 + 1)' ] ;                          
%                   
% Default input
f13dat.defval.Atr(natb).AttrName =  'sponge_generator_layer' ; 
f13dat.defval.Atr(natb).Unit = 'unitless' ;
f13dat.defval.Atr(natb).ValuesPerNode = 2 ;
f13dat.defval.Atr(natb).Val = [0.0 0] ;

if write == 1
    % write it out
    writefort13( f13dat, f13file )
end
%EOF
end

