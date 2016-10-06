function GRID_Z = DEM2GRD(DEM_X,DEM_Y,DEM_Z,GRID_X,GRID_Y,GRID_E,varargin)
% DEM2GRD: Uses the cell-averaged approach to interpolate the nodes
%          on the unstructured grid to the DEM. Ensure that inputs
%          are in projected coordinates (e.g. in metres)
%
%         Z = DEM2GRD(DEM_X,DEM_Y,DEM_Z,GRID_X,GRID_Y,'K',K,'SA',SA);
%         Input : DEM_X - The east-west matrix of the DEM grid of size
%                         m x n (m is number north-south points, n is
%                         number of east-west point)
%                 DEM_Y - The north-south matrix of the DEM grid of size
%                         m x n
%                 DEM_Z - Matrix of the depth or vertical difference from 
%                         the datum in the DEM of size m x n
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
%         SA (optional) - search radius. Helps to speed up the search 
%                         inside GRID_E by searching only adjacent node 
%                         numbers. Can get optimal results by renumbering 
%                         the mesh and using the node band width value as SA. 
%
%         Output : Z    - The depth or vertical difference from the datum
%                         in the mesh. length is nn
%
%	Author: William Pringle, CHL, Notre Dame University
%	Created: 2016-08-23
%
%   References:
%              V. Bilskie, S.C. Hagen (2013). "Topographic Accuracy 
%              Assessment of Bare Earth lidar-derived Unstructured Meshes".
%              Advances in Water Resources, 52, 165-177,
%              http://dx.doi.org/10.1016/j.advwatres.2012.09.003
%
%   Requires: bufferm2 - Function can be downloaded from:  
%   https://www.mathworks.com/matlabcentral/fileexchange/11095-bufferm2

%% Test optional arguments
if isempty(varargin)
    nn = length(GRID_X);
    K  = (1:nn)';        % K is all of the grid
    SA = int64(nn/10);   % Guess search radius
else
    if strcmp(varargin{1},'K')
        K  = varargin{2};
        nn = length(K);
    elseif strcmp(varargin{1},'SA')
        SA  = varargin{2};
    else
        disp(['argument ' varargin{1} ' invalid'])
        return
    end
    if length(varargin) == 4
        if strcmp(varargin{3},'K')
            K  = varargin{4};
            nn = length(K);
        elseif strcmp(varargin{3},'SA')
            SA  = varargin{4};
        else
            disp(['argument ' varargin{3} ' invalid'])
            return
        end
    else
        if strcmp(varargin{1},'K')
            SA  = int64(length(GRID_X)/10); % Guess search radius
        elseif strcmp(varargin{1},'SA')
            nn = length(GRID_X);
            K  = (1:nn)'; %K is all of the grid
        end
    end
end

%% Get grid area of DEM
DELTA_DEM = (DEM_X(1,2) - DEM_X(1,1)) ...
          * (DEM_Y(2,1) - DEM_Y(1,1));
      
%% Get average area of each node
% First obtain element areas
DELTA_E = polyarea(GRID_X(GRID_E(:,1:3))',GRID_Y(GRID_E(:,1:3))')';
ne = length(DELTA_E); DELTA_M = zeros(nn,1);
% Now average surrounding element areas for each node
for i = 1:nn
     if i == 1
        [I,~] = find(GRID_E(:,1:3) == K(i));
     else
        % search radius is SA
        ns = max(1,min(I)-SA);
        nee = min(ne,max(I)+SA);
        [I,~] = find(GRID_E(ns:nee,1:3) == K(i)); 
        I = I + double(ns) - 1;
     end
     if ~isempty(I)
        DELTA_M(i) = mean(DELTA_E(I));
     else
        % Search failed
        [I,~] = find(GRID_E(:,1:3) == K(i)); 
        DELTA_M(i) = mean(DELTA_E(I));
     end
end

%% Calculate CA - number of DEM grids to average - for each node
N   = 0.25*DELTA_M/DELTA_DEM; 
CA  = zeros(nn,1);
CA( N < 1)  = 1;
CA( N >= 1) = int64((2*N(N >= 1) + 1).^2);

%% Reshape DEM and delete uncessary parts
BufferL = sqrt(max(CA)*DELTA_DEM); % buffer of size sqrt(max(CA)*DELTA_DEM)
I = find(DEM_X(1,:) >= min(GRID_X(K)) - BufferL & ...
         DEM_X(1,:) <= max(GRID_X(K)) + BufferL);
J = find(DEM_Y(:,1) >= min(GRID_Y(K)) - BufferL & ...
         DEM_Y(:,1) <= max(GRID_Y(K)) + BufferL);
% Delete uncessary parts of DEM first step
DEM_X = DEM_X(J,I); DEM_Y = DEM_Y(J,I); DEM_Z = DEM_Z(J,I);  
DEM_s  = size(DEM_X); 
DEM_XX = reshape(DEM_X,DEM_s(2)*DEM_s(1),1); 
DEM_YY = reshape(DEM_Y,DEM_s(2)*DEM_s(1),1);
DEM_ZZ = reshape(DEM_Z,DEM_s(2)*DEM_s(1),1);
% Get rid of unecessary parts of DEM second step
CVH = convhull(GRID_X(K),GRID_Y(K)); % Get the convex hull
% Make buffer of size sqrt(max(CA)*DELTA_DEM)
[xb, yb] = bufferm2('xy',GRID_X(K(CVH)),GRID_Y(K(CVH)),BufferL,'out');
I = find(inpolygon(DEM_XX,DEM_YY,xb,yb) == 0);  % Check for DEM values outside convex hull with buffer
DEM_XX(I) = []; DEM_YY(I) = []; DEM_ZZ(I) = []; % Delete unecessary part of DEM

%% Get cell average for each node by averaging CA surrounding values in DEM
GRID_Z = zeros(nn,1);
% For only the nearest neighbour search
IDX = knnsearch([DEM_XX DEM_YY],[GRID_X(K(CA == 1)) GRID_Y(K(CA == 1))],'k',1);
GRID_Z(CA == 1) = DEM_ZZ(IDX);
% For cell-averaging searches
% Split search into 5 parts to avoid memory problems
kv = linspace(1,max(CA),5);
for k = 2:length(kv);
    kp = int64(kv(k-1));
    kn = int64(kv(k));
    % For less than or equal to kn and larger than kp
    IDX = knnsearch([DEM_XX DEM_YY],[GRID_X(K(CA > kp & CA <= kn)) ...
                                    GRID_Y(K(CA > kp & CA <= kn))],'k',kn); 
    if ~isempty(IDX)
        if size(IDX,1) == 1 || size(IDX,2) == 1
            GRID_Z(CA > kp & CA <= kn) = mean(DEM_ZZ(IDX));
        else
            GRID_Z(CA > kp & CA <= kn) = mean(DEM_ZZ(IDX),2)'; 
        end
    end
end
%EOF
end