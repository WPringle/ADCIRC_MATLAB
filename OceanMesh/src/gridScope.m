function [pc,tc,pinv,tinv]=gridScope(p,t)
tc = t; tinv = t ;
[~,bp]=extdom_edges2(t,p);

bbox = [min(bp(:,1)) max(bp(:,1))
    min(bp(:,2)) max(bp(:,2))];
bufx = 0.2*(bbox(1,2) - bbox(1,1));
bufy = 0.2*(bbox(2,2) - bbox(2,1));

figure; 
%simpplot(p,t);
% m_proj('Mercator','long',[bbox(1,1) - bufx, bbox(1,2) + bufx],...
%             'lat',[bbox(2,1) - bufy, bbox(2,2) + bufy])
% m_plot(bp(:,1),bp(:,2),'r.');
plot(bp(:,1),bp(:,2),'r.');
%m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
title('Scope out region...'); 


h    =  impoly;
poly =  h.getPosition;
poly = [poly; poly(1,:)];
hold on; plot(poly(:,1),poly(:,2),'k-','linewi',2);

bc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/+3; 
in =  inpoly(bc,poly);

tc(~in,:) = []; 
tidx = unique(tc(:)); 
%renumber tc to start at unity 
for i = 1 : length(tidx)
  g2l(tidx(i)) = i; 
end
for i = 1 : length(tc) 
   tc(i,1) = g2l(tc(i,1)); 
   tc(i,2) = g2l(tc(i,2)); 
   tc(i,3) = g2l(tc(i,3)); 
end
pc = p(tidx,:); 

g2l = []; 
tinv(in,:) = []; 
tidx = unique(tinv(:)); 
%renumber tinv to start at unity 
for i = 1 : length(tidx)
  g2l(tidx(i)) = i; 
end
for i = 1 : length(tinv) 
   tinv(i,1) = g2l(tinv(i,1)); 
   tinv(i,2) = g2l(tinv(i,2)); 
   tinv(i,3) = g2l(tinv(i,3)); 
end
pinv = p(tidx(:),:); 
clf, simpplot(pinv,tinv); 
end