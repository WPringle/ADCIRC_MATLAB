function B = SmoothBathSelection( EToV,xx,yy,B,searchr )
%SmoothBathSelection : Will plot a trimesh so you can select the desired 
%points using lasso. Then the function smooths these points using Gaussian
%distribution with specified search radius
%
% Inputs: EToV - n x 3 array of triangular elements (node indices)
%         xx   - n x 1 vector of x points (Cartesian)
%         yy   - n x 1 vector of y points (Cartesian)
%         B    - n x 1 vector of bathymetric depths
%      searchr - Radius of smoothing function (Cartesian)
%
% Outputs: B   - n x 1 vector of new bathymetric depths (with some
% smoothed)
%
% Useage:
% filename = 'fort.14';
% [EToV,VX,B,~,~,~] = readfort14(filename);
%
% CPP conversion:  
% lon0 = 75.214893 * pi/180; lat0 = -31.179511 * pi/180;
% R = 6378206.4;
% xx = VX(:,1) * pi/180;  yy = VX(:,2) * pi/180; 
% xx = R * (xx - lon0) * cos(lat0);
% yy = R * yy;
%
% B = SmoothBathSelection( EToV,xx,yy,B,1d3 );
%
% Requires: selectdata function at:
% https://www.mathworks.com/matlabcentral/fileexchange/13857-graphical-data-selection-tool
%
% BY William Pringle, Sep 15 2016

%% Define the guassian function
gaussC = @(dist, sigma2) (1 / sqrt(2 * pi )) * (exp(-dist./(2*sigma2)));   

%% Plot the mesh and ask for selection
if isempty(EToV)
   EToV = delaunay(xx,yy); 
end
trimesh(EToV,xx,yy,B)
view(2)
disp('Select points to smooth using lasso selection type')
[~,xs,ys] = selectdata('sel','Lasso'); 

%% Get the index of the selected points
K = zeros(size(xs));
for i = 1:length(xs)
   K(i) = find(xx == xs(i) & yy == ys(i)); 
end

%% smooth within search radius
[IDX, dist] = rangesearch([xx yy],[xs ys],searchr);
B1 = B(K);
for k = 1:length(K)
    % get gaussian smoothing function
    val = gaussC(dist{k}, 0.5*searchr);
    val = val/norm(val,1);
    B1(k) = val*B(IDX{k}); %sum the depths with weights
end
B(K) = B1;

end

