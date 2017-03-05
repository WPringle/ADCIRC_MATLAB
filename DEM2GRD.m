function GRID_Z = DEM2GRD(DEM_X,DEM_Y,DEM_Z,GRID_X,GRID_Y,GRID_E,varargin)
% DEM2GRD: Uses the cell-averaged approach to interpolate the nodes
%          on the unstructured grid to the DEM. Ensure that inputs
%          are in projected coordinates (e.g. in metres)
%
%         Z = DEM2GRD(DEM_X,DEM_Y,DEM_Z,GRID_X,GRID_Y,'K',K,'SA',SA);
%         Input : DEM_X - The east-west matrix of the DEM grid
%                         in ndgrid format
%                 DEM_Y - The north-south matrix of the DEM grid
%                         in ndgrid format
%                 DEM_Z - Matrix of the depth or vertical difference from 
%                         the datum in the DEM in ndgrid format
%         
%                GRID_X - A vector of the east-west locations of all the
%                         nodes in the mesh, length in number of nodes, nn
%                GRID_Y - A vector of the north-south locations of all the
%                         nodes in the mesh, length is nn
%                GRID_E - Matrix size ne x 3 of the node distribution of
%                         all the elements, ne of the mesh
%
%          K (optional) - vector of relevant nodes to search. This can
%                         significantly speed up the calculation by either 
%                         only interpolating part of the mesh or 
%                         intepolating whole mesh but in segments. Function
%                         automatically removes uncessary portions of DEM
%                         Example to call: 
%                         K = find( GRID_X >= lon_min & GRID_X <= lon_max);
%                         Z(K) = DEM2GRD(DEM_X,DEM_Y,DEM_Z,GRID_X,GRID_Y,'K',K);
%
%         Output : Z    - The depth or vertical difference from the datum
%                         in the mesh. length is nn
%
%	Author: William Pringle, CHL, Notre Dame University
%	Created: 2016-08-23
%   Updated: 2017-02-14, DELTA_M algorithm improved without searching reqd, 
%                        correction of CA averaging at end.
%   Updated: 2017-03-05, fixed CA formula, needed a sqrt.
%
%   References:
%              V. Bilskie, S.C. Hagen (2013). "Topographic Accuracy 
%              Assessment of Bare Earth lidar-derived Unstructured Meshes".
%              Advances in Water Resources, 52, 165-177,
%              http://dx.doi.org/10.1016/j.advwatres.2012.09.003
%
%   Requires: bufferm2 - Function can be downloaded from:  
%   https://www.mathworks.com/matlabcentral/fileexchange/11095-bufferm2
%   *This is actually now commented out in line 93 and thus is not necessarily reqd.

%% Test optional arguments
if isempty(varargin)
    nn = length(GRID_X);
    K  = (1:nn)';        % K is all of the grid
elseif strcmp(varargin{1},'K')
    K  = varargin{2};
    nn = length(K);
else
    disp(['argument ' varargin{1} ' invalid'])
    return
end

%% Get grid area of DEM
DELTA_DEM = (DEM_X(2,1) - DEM_X(1,1)) ...
          * (DEM_Y(1,2) - DEM_Y(1,1));
      
%% Get average area of each node
% First obtain element areas
disp('Getting average element area for each node to calculate CA')
DELTA_E = polyarea(GRID_X(GRID_E(:,1:3))',GRID_Y(GRID_E(:,1:3))')';
ne = length(DELTA_E); nnn = length(GRID_X); 
DELTA_M = zeros(nnn,1); COUNT = zeros(nnn,1);
%% Sum area over all nodes of that element and divide by count
for i = 1:ne
    DELTA_M(GRID_E(i,:)) = DELTA_M(GRID_E(i,:)) + DELTA_E(i);
    COUNT(GRID_E(i,:)) = COUNT(GRID_E(i,:)) + 1;
end
DELTA_M = DELTA_M(K)./COUNT(K);

%% Calculate CA - number of DEM grids to average - for each node
N   = 0.25*sqrt(DELTA_M/DELTA_DEM); 
CA  = zeros(nn,1);
CA( N < 1)  = 1;
CA( N >= 1) = int64((2*N(N >= 1) + 1).^2);
%% Delete uncessary parts of DEM
BufferL = sqrt(max(CA)*DELTA_DEM); % buffer of size sqrt(max(CA)*DELTA_DEM)
I = find(DEM_X(:,1) >= min(GRID_X(K)) - BufferL & ...
         DEM_X(:,1) <= max(GRID_X(K)) + BufferL);
J = find(DEM_Y(1,:) >= min(GRID_Y(K)) - BufferL & ...
         DEM_Y(1,:) <= max(GRID_Y(K)) + BufferL);
% Delete uncessary parts of DEM first step
DEM_X = DEM_X(I,J); DEM_Y = DEM_Y(I,J); DEM_Z = DEM_Z(I,J);  
% % Get rid of unecessary parts of DEM second step
% CVH = convhull(GRID_X(K),GRID_Y(K)); % Get the convex hull
% % Make buffer of size sqrt(max(CA)*DELTA_DEM)
% [xb, yb] = bufferm2('xy',GRID_X(K(CVH)),GRID_Y(K(CVH)),double(BufferL),'out');
% I = find(InPolygon(DEM_X,DEM_Y,xb,yb) == 0);  % Check for DEM values outside convex hull with buffer
% DEM_X(I) = []; DEM_Y(I) = []; DEM_Z(I) = []; % Delete unecessary part of DEM

%% Get cell average for each node by averaging CA surrounding values in DEM
GRID_Z = zeros(nn,1); 
% Get gridded interpolant where CA == 1
% For only the gridded Interpolant interpolation
disp('doing spline interpolation where CA = 1')
F = griddedInterpolant(DEM_X,DEM_Y,DEM_Z,'spline');
GRID_Z(CA == 1) = F(GRID_X(K(CA == 1)),GRID_Y(K(CA == 1)));

% Return if not CA is larger than 1
if isempty(CA(CA>1))
    return;
end
% For cell-averaging searches where CA > 1   
% Reshape for kD-tree
DEM_s  = size(DEM_X); 
DEM_X = reshape(DEM_X,DEM_s(2)*DEM_s(1),1); 
DEM_Y = reshape(DEM_Y,DEM_s(2)*DEM_s(1),1);
DEM_Z = reshape(DEM_Z,DEM_s(2)*DEM_s(1),1);

% set up searching algorithm
disp('setting up knnsearcher for CA > 1')
Mkd = KDTreeSearcher([DEM_X DEM_Y],...
                      'BucketSize',ceil(100^(2/3)*median(CA(CA>1))^(1/3)));
CA_un = unique(CA,'sorted'); CA_un(1) = [];
disp('doing the search for each value of CA')
for k = CA_un'
    % Use kd-tree
    IDX = knnsearch(Mkd,[GRID_X(K(CA == k)) GRID_Y(K(CA == k))],'k',k);
    GRID_Z(CA == k) = mean(DEM_Z(IDX),min(2,size(IDX,1))); 
end
%EOF
end