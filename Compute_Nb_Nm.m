function [Nb,Nm] = Compute_Nb_Nm(lon_M,lat_M,B,zcontour,...
                                 N,lon_N,lat_N,lon0,lat0)
% Compute_Nb_Nm: Compute the buoyancy frequency, N at the seabed (Nb) and 
%                the mean over the depth (Nm) for an unstructured mesh from
%                scattered grid points of N at specified contours
%
% [Nb,Nm] = Compute_Nb_Nm(lon_M,lat_M,B,zcontour,N,lon_N,lat_N,lon0,lat0)
% Input : lon_M    - longitude points of nodes in mesh
%         lon_M    - latitude points of nodes in mesh
%         B        - depths of nodes in mesh
%         zcontour - the contours where we have values of N
%         N        - cell of buoyancy freq. scatter points for each contour
%         lon_N    - cell of longitude scatter points for each contour
%         lat_N    - cells of latitude scatter points for each contour
%         lon0     - longitude centre of CPP conversion
%         lat0     - latitude centre of CPP conversion 
%
% Output : Nb     - Buoyancy frequency at seabed
%          Nm     - Depth-averaged buoyancy frequency over depth
%
% Author: William Pringle, CHL, Notre Dame University
% Created: 2016-10-05

% Set minimum length for interpolating data
Min_length = 30;

% Get xx and yy by CPP conversion
R = 6378206.4; % radius of earth

xx = lon_M * pi/180;  yy = lat_M * pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

%% Calculation
% initialisation
Nb = zeros(size(B)); Nm = zeros(size(B)); 
Nlast = zeros(size(B)); Count = zeros(size(B));
contour_int = zcontour(2) - zcontour(1);
% loop
for zvalue = 1:length(zcontour)

    if length(lat_N{zvalue}) < Min_length; continue; end
    
    x = lon_N{zvalue};
    y = lat_N{zvalue};
    x = double(x); y = double(y);
    % CPP conversion
    x = x * pi/180;  y = y * pi/180; 
    x = R * (x - lon0) * cos(lat0);
    y = R * y; 
    
    % Make the interpolant using scatteredInterpolant natural
    F = scatteredInterpolant(x,y,N{zvalue},'natural','nearest');
    
    % Get data where depth is equal to rounded zvalue
    I = find( contour_int*int32(B/contour_int) == zcontour(zvalue));
    
    % Interpolate data to get Nb 
    Nb(I) = F(xx(I),yy(I));
    
    % Get data where depth is larger than or equal to round zvalue
    J = find( contour_int*int32(B/contour_int) >= zcontour(zvalue));
    
    % Add data at each contour for averaging later
    Nlast(J) = F(xx(J),yy(J));
    Nm(J) = Nm(J) + Nlast(J);
    Count(J) = Count(J) + 1;
end
% Divide by the count of contours
Nm(Count > 0) = Nm(Count > 0)./Count(Count > 0);
% Constant extrapolation to depths outside of last contour 
Nb(Nb == 0 & Count > 0) = Nlast(Nb == 0 & Count > 0);
%EOF
end