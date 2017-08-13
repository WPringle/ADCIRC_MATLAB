function [adj,xadj]=EleToEle(t)
%--if an element shares an edge then it is a neighbor 
%--returns element-to-element connectivity in crs format. 
%--triangle ie is connected to adj(xadj(ie):xadj(ie+1)-1,1) triangles
% kjr 2017, fixed August 2017 for self adj triangles (i.e., disjoint
% elements). 
t = sort(t,2); 
edges = [t(:,[1 2]);t(:,[1 3]);t(:,[2 3])];
nt = size(t,1);
trinum = repmat((1:nt)',3,1);
[~,tags] = sort(edges*[2^31;1]);
edges=edges(tags,:); 

trinum = trinum(tags);
k = find(all(diff(edges,1)==0,2));
adj=trinum([k,k+1]);

[~,dmy1]=sort(adj(:,1));
[~,dmy2]=sort(adj(:,2)); 

junk=[adj(dmy1,:); fliplr(adj(dmy2,:)); (1:nt)',(1:nt)'];%kjr added self adjs. 
[~,idx2]=sort(junk*[2^31;1]);  
adj=junk(idx2,:); 

xadj = find(diff(adj(:,1))==1); 
adj(:,1)=[];
xadj=[xadj;length(adj)]; 
xadj = [1;xadj+1];                                              
end