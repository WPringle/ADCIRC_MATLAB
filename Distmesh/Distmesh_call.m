% Distmesh calling script
clearvars; close all; 

%% Setting parth and compling InPolygon if required
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB')
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB/distmesh1.2')
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB/m_map')
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB/InPolygon-MEX')
%mex '/afs/crc.nd.edu/user/w/wpringle/MATLAB/InPolygon-MEX/InPolygon.c' 

%% Set up parallel pool
Proc_num = 4;
if isempty(gcp)
    parpool('local',Proc_num);
end

%% Set parameters
mapfile = 'IDIOMS_v7_split.map';                           
bathyfile = ['E:\Global_Data\SRTM30_PLUS_w_Abyssal_Hills\' ...
            'bathy_SSG_1_120_GLOBAL_landmask.nc'];
edgelength = 1/120;    % min edgelength in degrees
dist_param = 0.1;      % Distance paramater (percent that the edgelength 
                       % should change with distance, set zero to ignore)
wl_param = 240;        % parameter in wavelength function (set zero to ignore)
slope_param = 30;      % parameter in slope function (set zero to ignore)
itmax       = 1000;    % Maximum number of iterations allowed in distmesh
plot_on = 1;           % Plot? (Yes = 1, No = 0)
%finame = 'Distmesh_2km/split_2km_imp';
outfiname  = 'distmesh_1km_split/split_1km'; % Output .mat name;

% Load split up mesh
%load([finame '_' num2str(split) '.mat'])
ini_p = [];
pfix = [];
if exist('p','var')
    ini_p = p;
end
% [EToV,p,~,opedat,boudat,title] = readfort14( finame );
% %Note about finding edges - must not have any nodes that are used in more
% %than two edges - i.e. apart of more than one polygon each
% [etbv,vxe,~,~] = extdom_edges( EToV, p ) ;
% 
% % Give the end points for gridscope boundary
% fixp_s = [56.8144047817, 24.2723540496]; 
% fixp_s1 = [56.8694017838, 24.3596682147]; 
% fixp_e = [56.8675699926, 27.0544989123];
% kdx = knnsearch(vxe,[fixp_s;fixp_s1;fixp_e]);
% % Find a segment that is made up of [ik1 ik2]
% idx = unique(etbv(:)) ; 
% iseg = sort([idx(kdx(1)) ; idx(kdx(2))])*ones( 1, length(idx)) ;
% iedbeg = find(~(sum(sort(etbv, 1) - iseg))) ;
% for ipsbeg = 1:2
%     if ipsbeg == 1
%         ipsend = 2; %
%     else
%         ipsend = 1; %
%     end
%     [vso,idv,ide] = extdom_polygon( etbv, p, iedbeg, ipsbeg, ipsend ) ;
%     % vso(:,1:2) - coordinate of an extractred polygon 
%     % idv  -- indices of extracted ploygon in the global mesh,  vso = pv(idv,:)  
%     % ide  -- indeces of edges in etbv that constitute vso,  ie. vso = vxe(reshape(etbv(:,ide),length(ide)*2),:)  
%     cw = ispolycw(vso(:,1), vso(:,2));
%     if ~cw; break; end
% end
% 
% kdx = knnsearch(p,[fixp_s;fixp_s1;fixp_e]);
% idvs = find(idv == kdx(1));
% idve = find(idv == kdx(3));
% % Find the gridscope boundary
% fixp = vso(idvs(1):idve(1),:);

%% Call distmesh
for split = 1:26
    %if split == 12 || split == 17; continue; end
    tic
    General_distmesh(mapfile,bathyfile,edgelength,dist_param,...
                          wl_param,slope_param,itmax,plot_on,ini_p,pfix,split);
    toc
    disp(['finished split no.' num2str(split)])
    %% Save outputs to .mat file
    %writefort14( outfiname, t, p, zeros(length(p),1), opedat , boudat,title)
%     if ~isempty(p)
%         save([outfiname '_' num2str(split) '.mat'],'p','t');
%     end
end
%% Close parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

