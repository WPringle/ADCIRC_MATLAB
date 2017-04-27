%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete not connected elements                                           %
% This script does a spider search through all elements and gets rid      %
% of the elements not connected to the main mesh                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all;

% Load mesh (enter mesh name with suffix)
fort14name = 'South_Carolina3';
[t,p,~,~,~,~] = readfort14([fort14name '.grd']);

% Give element on boundary open boundary
EToS =  508905;

% Get connectivity
EToE = Connect2D(t);

%% Traverse grid deleting elements outside
ic = zeros(ceil(sqrt(length(t))*2),1);
ic0 = zeros(ceil(sqrt(length(t))*2),1);
nflag = zeros(length(t),1);
icc0 = 1;
ic0(1) = EToS;
% loop until convergence is read
while 1
    icc = 0;
    for nn = 1:icc0
        i = ic0(nn);
        % Flag the current element as OK
        nflag(i) = 1;
        % Search neighbouring elements
        for nb = EToE(i,:)
           if nflag(nb) == 0
               icc = icc + 1;
               ic(icc) = nb;
               nflag(nb) = 1;
           end
        end
    end
    if icc ~= 0 
        icc0 = icc;
        ic0(1:icc0) = ic(1:icc0);
    else
        break
    end
end

% Only keep nflag == 1
t = t(nflag == 1,:);
% delete disjoint nodes
[p,t] = fixmesh(p,t);

%% Mesh over small polygons inside mesh
% Get boundaries
%[etbv,vxe,~,~] = extdom_edges( t, p ) ;
%
%iedbeg = 1;
%[vso,idv,ide,error] = extdom_polygon( etbv, p, iedbeg, 1, 2 ) ;

% plot
simpplot(p,t)
