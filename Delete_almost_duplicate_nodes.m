% Find almost duplicate nodes, delete those and rearrange elements
% Sometimes SMS creates duplicate nodes when merging meshes but they do not
% appear as actual duplicate nodes as they are at slightly different
% location

clearvars; clc;

% set small distance
d_small = 1e-5;

% load the mesh
load IDIOMS_v6.3.mat

% Knnsearch
[IDX,DIST] = knnsearch(VX,VX,'k',2);
% delete own node
IDX(:,1) = [];
DIST(:,1) = [];

% Find very small distances
IA = find(DIST < d_small);

% Loop through and replace all relevant indices in EToV
n = 0;
while n < length(IA)
   n = n + 1;
   node_n = IA(n);
   neigh_n = IDX(IA(n));
   [I,J] = find(EToV == neigh_n); 
   for ii = 1:length(I)
        EToV(I(ii),J(ii)) = node_n;
   end
   IA( IA == neigh_n ) = [];
end

% Find large enough distances
IC = find(DIST >= d_small);

% Add the remaining IA to IC and sort
IA = [IA; IC];
IA = sortrows(IA);

% write out with numbering preserved
writefort14_preservenumbering( 'IDIOMS_v6.3_correct.grd', EToV, IA, ...
                               VX(IA,:), B(IA), opedat, boudat, title)

