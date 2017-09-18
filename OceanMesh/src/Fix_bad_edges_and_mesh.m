function [ p , t ] = Fix_bad_edges_and_mesh( p , t, nscreen )
%  This script first checks and deletes for nodes that are connected to
%  more than 2 boundary edges so that the mesh is well-formed. The it finds
%  small disjoint portions of the graph and removes them using a
%  "spider-search". The disjoint portions are kept given that their
%  composition represents less than 5% of the total graph. The program
%  alternates between checking interior and exterior portions of the graph
%  exhaustively. 
%  Written by William Pringle, 
%    modifications by Keith Roberts, CHL,UND 2017

% Ensure mesh is "fixed"
[p,t] = fixmesh(p,t);

[p,t] = delete_elements_inside_main_mesh(p,t,nscreen);

disp('ALERT: finished cleaning up mesh..'); 

end

function [p,t] = delete_elements_outside_main_mesh(p,t,nscreen)
%% Delete all elements outside the main mesh
t1 = t; t = []; L = length(t1); close all;
max_lim = 10; % sometimes the random element was a bad choice in this case, we may never converge. break out and try again in this case. kjr
while 1
    % Give a random element not outside main mesh
    EToS =  randi(length(t1),1);
    min_del = L;
    counter = 0;
    while 1
        % Get connectivity if it changed or first time.
        [EToE,xadj] = EleToEle(t1);
        
        % Traverse grid deleting elements outside
        ic = zeros(ceil(sqrt(length(t1))*2),1);
        ic0 = zeros(ceil(sqrt(length(t1))*2),1);
        nflag = zeros(length(t1),1);
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
                for nb = EToE(xadj(i):xadj(i+1)-1,2)'
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
        if nscreen
            disp(['  CANDIDATE: deleting ' num2str(length(find(nflag == 0))) ...
                ' elements outside main mesh'])
        end
        
        counter = counter + 1;
        if(counter >= max_lim)
            disp('ALERT: TOO MANY ITERATIONS...RESTARTING');
            break; 
        end; 
        
        if nnz(find(nflag == 0))/length(t1) < 0.5 || ...
                nnz(find(nflag == 0)) == min_del
            % choice of EToS OK
            % (deleting less than half of the triangulation)
            break
            
        else
            min_del = min(min_del,nnz(find(nflag == 0)));
            EToS = find(nflag == 0); EToS = EToS(1);
        end
    end
    % adding to the triangulation
    t = [t; t1(nflag == 1,:)];
    % deciding whether portion is small enough to exit loop or not
    if nnz(find(nflag == 0))/L < 0.05
        disp(['ACCEPTED: deleting ' num2str(length(find(nflag == 0))) ...
            ' elements outside main mesh']) ;
        break
    else
        if(nscreen)
            disp('  REJECTED CANDIDATE..retrying');
        end
        % making the subset triangulation
        t1 = t1(nflag == 0,:);
    end
end

% delete disjoint nodes
[p,t] = fixmesh(p,t);

end

function [p,t] = delete_elements_inside_main_mesh(p,t,nscreen)
%% Delete some elements so that we have no nodes
%%  connected to more than 2 boundary edges
% loop until no longer exist
while 1
    % Get boundaries
    [etbv,vxe] = extdom_edges2( t, p ) ;
    if length(etbv) == length(vxe); break; end
    %
    % Get all nodes that are on edges
    [nodes_on_edge,~,n] = unique (etbv(:));
    % Count how many edges a node appears in
    I = accumarray(n,1:numel(n),[],@(x){x});
    count=cellfun('length',I);
    %
    [vtoe,nne] = VertToEle(t);
    % Get the nodes which appear more than twice and delete element connected
    % to this nodes where all nodes of element are on boundary edges
    del_elem_idx = [];
    for i = nodes_on_edge(count > 2)'
        con_elem = vtoe(1:nne(i),i);
        n = 0; del_elem = [];
        for elem = con_elem'
            %I = find(etbv(:) == t(elem,1), 1);
            %J = find(etbv(:) == t(elem,2), 1);
            %K = find(etbv(:) == t(elem,3), 1);
             I = etbv(:) == t(elem,1); 
             J = etbv(:) == t(elem,2); 
             K = etbv(:) == t(elem,3); 

            % all nodes on element are boundary edges
            %if ~isempty(I) && ~isempty(J) && ~isempty(K)
            if any(I) && any(J) && any(K) 
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
        disp(['ACCEPTED: deleting ' num2str(length(del_elem_idx)) ...
            ' elements inside main mesh'])
    t(del_elem_idx,:) = [];
    
    % delete disjoint nodes
    [p,t] = fixmesh(p,t);
    
    % Delete elements outside to ensure covergence
    [p,t] = delete_elements_outside_main_mesh(p,t,nscreen);
end

end


