%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fix bad mesh elements
% This script loops through the boundaries of the mesh and fixes up
% disconnected boundaries in the following ways. 
% 1) any elements outside of the major outer polygon are removed
% 2) non-elements inside the major outer polygon are turned into elements
% 3) any elements inside of the island polygons are removed
% 4) non-elements outside the island polygon are turned into elements
%
% At times you may be asked to make a decision on which route to traverse 
% along the boundary. Enter 1, 2 or 3 as prompted. The possible routes to
% choose from will be shown along with the surrounding node and element
% setup so that you can make a decision. Such decisions arise when the
% algorithm cannot easily eliminate enough bad routes.  

clearvars; close all;

% Load mesh (enter mesh name with suffix)
fort14name = 'WESPAC_112_5';

[ev,pv,~,~,~,~] = readfort14([fort14name '.grd']);

% Give a random node number of the major outer polygon
ob_se = 1;

% Extrace edges on the boundaries from a given mesh, pv, ev %
[etbv,vxe,etoe,etof] = extdom_edges( ev, pv ) ;

if exist('ob_se','var')
    iedbeg  = knnsearch(pv(etbv(1,:),:),pv(ob_se,:));
else
    iedbeg  = 1 ;
end
start = 1;
ipsbeg = 1 ; % polygon travese in 'ipsbeg' --> 'ipsend' 
ipsend = 2 ; %
%% Extract the outter domain 
polygon.outer = []; polygon.inner = [];
while ~isempty(etbv)
    if start == 0
       iedbeg = 1; 
    end  
    [vso,idv,ide] = extdom_polygon_fix( etbv, ev, pv, iedbeg, ipsbeg, ipsend ) ;
    etbv(:,ide) = []; 
    if length(vso) < 7; continue; end
    figure(2);
    plot(vso(:,1),vso(:,2))
    hold on
    if start == 1
        polygon.outer = vso;
        start = 0;
        % Get stuff outside of the current polygon
        In = InPolygon(pv(etbv(1,:),1),pv(etbv(1,:),2),vso(:,1),vso(:,2));
        etbv(:,In == 0) = [];
        In = InPolygon(pv(etbv(2,:),1),pv(etbv(2,:),2),vso(:,1),vso(:,2));
        etbv(:,In == 0) = [];
    else
        if vso(1,1) ~= vso(end,1)
           vso(end+1,:) = vso(1,:); 
        end
        cw = ispolycw(vso(:,1),vso(:,2));
        if cw 
            % we are clockwise, make it acw
            polygon.inner = [polygon.inner; NaN NaN; flipud(vso)]; 
        else
            % we are anticlockwise
            polygon.inner = [polygon.inner; NaN NaN; vso]; 
        end
        % Get stuff inside of the current polygon
        In = InPolygon(pv(etbv(1,:),1),pv(etbv(1,:),2),vso(:,1),vso(:,2));
        etbv(:,In == 1) = [];
        In = InPolygon(pv(etbv(2,:),1),pv(etbv(2,:),2),vso(:,1),vso(:,2));
        etbv(:,In == 1) = [];
    end
end
%%
% Now we have all polygons lets get the inpolygon and delete all elements
% outside (for outer) and inside (for inner) of this
pmid = (pv(ev(:,1),:) + pv(ev(:,2),:) + pv(ev(:,3),:))/3;
In = InPolygon(pmid(:,1),pmid(:,2),polygon.outer(:,1),polygon.outer(:,2));
ev(In == 0,:) = [];
pmid = (pv(ev(:,1),:) + pv(ev(:,2),:) + pv(ev(:,3),:))/3;
In = InPolygon(pmid(:,1),pmid(:,2),polygon.inner(:,1),polygon.inner(:,2));
ev(In == 1,:) = [];
%
IA = unique(ev);
pv = pv(IA,:);
%
% Now we have new points lets do the delaunay to make triangles inside the
% polygon that we don't currently have
ev = delaunay(pv(:,1),pv(:,2));
% delete shitty triangles again..
pmid = (pv(ev(:,1),:) + pv(ev(:,2),:) + pv(ev(:,3),:))/3;
In = InPolygon(pmid(:,1),pmid(:,2),polygon.outer(:,1),polygon.outer(:,2));
ev(In == 0,:) = [];
pmid = (pv(ev(:,1),:) + pv(ev(:,2),:) + pv(ev(:,3),:))/3;
In = InPolygon(pmid(:,1),pmid(:,2),polygon.inner(:,1),polygon.inner(:,2));
ev(In == 1,:) = [];
%
% Finalise mesh - deleting disjoint nodes
[pv,ev] = fixmesh(pv,ev);
% write 
writefort14( [fort14name '_imp.grd'] , ...
              ev, pv, zeros(length(pv),1), [] , [] ,'grid' ) ;