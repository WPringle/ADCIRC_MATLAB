function [etoe,idx]=EleToEle(t)
%--if an element shares an edge then it is a neighbor 
%--returns element-to-element connectivity in crs format. 
%--triangle ie is connected to etoe(idx(ie):idx(ie+1)-1,2) triangles
% kjr 2017
t = sort(t,2);
edges = [t(:,[1 2]);t(:,[1 3]);t(:,[2 3])];
nt = size(t,1);
trinum = repmat((1:nt)',3,1);
[edges,tags] = sortrows(edges);
trinum = trinum(tags);
k = find(all(diff(edges,1)==0,2));
etoe=trinum([k,k+1]);
etoe=sortrows(etoe);                         
idx = find(diff(etoe)==1); idx = [1;idx+1];                                              
end