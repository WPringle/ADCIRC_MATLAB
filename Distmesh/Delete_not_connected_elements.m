%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delete not connected elements                                           %
% This script does a spider search through all elements and gets rid      %
% of the elements not connected to the main mesh                          %
% Then it will search for bad connected edges inside main mesh and delete %
% elements strategically until no more bad connected edges exist          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all;

% Load mesh (enter mesh name with suffix)
% fort14name = 'South_Carolina3';
% [t,p,~,~,~,~] = readfort14([fort14name '.grd']);
load New_WesPac_Ocean_0_smooth

% Give a random element not outside main mesh
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

disp(['deleting ' num2str(length(find(nflag == 0))) ...
      ' elements outside main mesh'])
% Only keep nflag == 1
t = t(nflag == 1,:);
% delete disjoint nodes
[p,t] = fixmesh(p,t,1d-6);

simpplot(p,t)
%% Delete some elements so that we have no nodes 
%%  connected to more than 2 boundary edges
% loop until no longer exist
while 1
    % Get boundaries
    [etbv,vxe,~,~] = extdom_edges( t, p ) ;
    if length(etbv) == length(vxe); break; end
    %
    % Get all nodes that are on edges
    nodes_on_edge = unique(etbv(:));
    %
    % Count how many edges a node appears in
    N = numel(nodes_on_edge);
    count = zeros(N,1);
    for k = 1:N
        count(k) = sum(etbv(:)==nodes_on_edge(k));
    end
    %
    [ vtov, nnv, vtoe, nne, ~ ] = NodeConnect2D( t );
    % Get the nodes which appear more than twice and delete element connected
    % to this nodes where all nodes of element are on boundary edges
    del_elem_idx = [];
    for i = nodes_on_edge(count > 2)'
        con_elem = vtoe(1:nne(i),i);
        n = 0; del_elem = [];
        for elem = con_elem'
           I = find(etbv(:) == t(elem,1));
           J = find(etbv(:) == t(elem,2));
           K = find(etbv(:) == t(elem,3));
           % all nodes on element are boundary edges
           if ~isempty(I) && ~isempty(J) && ~isempty(K)
               n = n + 1;
               del_elem(n) = elem;
           end
        end
        if n == 1
           del_elem_idx(end+1) = del_elem;
        elseif n > 1
           tq = gettrimeshquan( p, t(del_elem,:));
           [~,idx] = min(tq.qm);
           % delete shittiest element
           del_elem_idx(end+1) = del_elem(idx);
        else
           % no connected elements have all nodes on boundary edge so we
           % just pick a random shitty element to delete
           tq = gettrimeshquan( p, t(con_elem,:));
           [~,idx] = min(tq.qm);
           % delete shittiest element
           del_elem_idx(end+1) = con_elem(idx);
           %del_elem_idx(end+1:end+length(con_elem)) = con_elem';
        end
    end
    disp(['deleting ' num2str(length(del_elem_idx)) ...
          ' elements inside main mesh'])
    t(del_elem_idx,:) = [];
    % delete disjoint nodes
    [p,t] = fixmesh(p,t);
end
disp('finished cleaning up mesh, plotting and saving..')
% plot
simpplot(p,t)
%
save('New_WesPac_Ocean_0_smooth_imp.mat','p','t')