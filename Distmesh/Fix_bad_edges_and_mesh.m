function [ p , t ] = Fix_bad_edges_and_mesh( p , t )
% [ p , t ] = Fix_bad_edges_and_mesh( p , t )
%   Detailed explanation goes here

% Ensure mesh is "fixed"
[p,t] = fixmesh(p,t);

% First delete unecessary elements outside main mesh (fast)
[p,t] = delete_elements_outside_main_mesh(p,t);

% Second, delete unecessary elements inside main mesh (iterative; slow)
[p,t] = delete_elements_inside_main_mesh(p,t);

% Delete uncessary elements outside main mesh again (fast)
[p,t] = delete_elements_outside_main_mesh(p,t);

disp('finished cleaning up mesh..')

% Now move nodes in bad remaining elements (fast)
[p,t] = fix_interior_angles(p,t);

disp('finished moving nodes in bad elements..')
end

function [p,t] = delete_elements_outside_main_mesh(p,t)
%% Delete all elements outside the main mesh
while 1
    % Give a random element not outside main mesh
    EToS =  randi(length(t),1);

    % Get connectivity
    EToE = Connect2D(t);

    % Traverse grid deleting elements outside
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

    if length(find(nflag == 0))/length(t) < 0.5
       % choice of EToS OK
       break
    end
end
  
% Only keep nflag == 1
t = t(nflag == 1,:);

% delete disjoint nodes
[p,t] = fixmesh(p,t); %,1d-6);

end

function [p,t] = delete_elements_inside_main_mesh(p,t)
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
    [ ~, ~, vtoe, nne, ~ ] = NodeConnect2D( t );
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

end

function [p,t] = fix_interior_angles(p,t)
%% Identify bad elements and move one node to centre of surrounding ones
    [ vtov, nnv, ~, ~, ~ ] = NodeConnect2D( t );
    tq = gettrimeshquan( p, t);
    [I,~] = find(tq.vang < 30*pi/180 | tq.vang > 130*pi/180);
    %       
    for i = I'
        % Get the node belonging to largest angle
        [ang,j] = max(tq.vang(i,:));
        if ang < pi*90/180
            % Get node belonging to smallest angle
            [~,j] = min(tq.vang(i,:));
        end
        node = t(i,j) ;
        % Move this node to center of surrounding nodes
        center = mean(p(vtov(1:nnv(node),node),:));
        p(node,:) = center;
    end
end