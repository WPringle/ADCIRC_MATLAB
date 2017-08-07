function  [bnde,bnd] = extdom_edges2( t , p )
%
% Given Element table, extract boundary edges
edge = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];                               % Non-unique edges
edge = sortrows(sort(edge,2));                                             % Put shared edges next to each other
idx = all(diff(edge)==0,2);                                                % Find shared edges
idx = [idx;false]|[false;idx];                                             % True for all shared edges
bnde = edge(~idx,:);                                                       % Boundary edges
%edge = edge(idx,:);                                                       % Internal edges
%edge = [bnde; edge(1:2:end-1,:)];                                         % Unique edges
bnd  = unique(bnde);                                                       % Boundary nodes
bnde = bnde; bnd=p(bnd,:); 
end
