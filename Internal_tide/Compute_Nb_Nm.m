function [Nb,Nm] = Compute_Nb_Nm(lon_M,lat_M,B,zcontour,N,lon_N,lat_N,proj)
                                  %,lon0,lat0)
% Compute_Nb_Nm: Compute the buoyancy frequency, N at the seabed (Nb) and 
%                the mean over the depth (Nm) for an unstrcutred mesh from
%                scattered grid points of N at specified contours
%
% [Nb,Nm] = Compute_Nb_Nm(lon_M,lat_M,B,zcontour,N,lon_N,lat_N,lon0,lat0)
% Input : lon_M    - longitude points of nodes in mesh
%         lat_M    - latitude points of nodes in mesh
%         B        - depths of nodes in mesh
%         zcontour - the contours where we have values of N
%         N        - cell of buoyancy freq. scatter points for each contour
%         lon_N    - cell of longitude scatter points for each contour
%         lat_N    - cells of latitude scatter points for each contour
%         proj     - string that defines the projection to use
%                    e.g. 'Mercator'.., (type m_proj('set') for options)
%
% Output : Nb     - Buoyancy frequency at seabed
%          Nm     - Depth-averaged buoyancy frequency over depth
%
% Author: William Pringle, CHL, Notre Dame University
% Created: 2016-10-05

% Set minimum length for interpolating data
Min_length = 20;

% Projection set
% Getting minima and maxima
loncmin = cellfun(@min, lon_N,'UniformOutput',0);
loncmax = cellfun(@max, lon_N,'UniformOutput',0);
latcmin = cellfun(@min, lat_N,'UniformOutput',0);
latcmax = cellfun(@max, lat_N,'UniformOutput',0);
loncmin = min([loncmin{:}]);
loncmax = max([loncmax{:}]);
latcmin = min([latcmin{:}]);
latcmax = max([latcmax{:}]);
% Ignore data outside of a large enough domain 
lon_min = double(max(loncmin,min(lon_M(:))-30));
lon_max = double(min(loncmax,max(lon_M(:))+30));
lat_min = double(max(latcmin,min(lat_M(:))-30));
lat_max = double(min(latcmax,max(lat_M(:))+30));
if lat_max == 90; lat_max = 89; end
m_proj(proj,'lon',[ lon_min lon_max],'lat',[ lat_min lat_max])

% Conversion to projection coordinates
[xx,yy] = m_ll2xy(lon_M,lat_M);        
% 
xx = double(xx);
yy = double(yy);

%% Calculation
% initialisation
Nb = zeros(size(B)); Nm = zeros(size(B)); 
Nlast = zeros(size(B)); Count = zeros(size(B));
% do the interpolation onto the mesh
N_interp = cell(length(zcontour),1);
for zvalue = 1:length(zcontour)

    x = lon_N{zvalue};
    y = lat_N{zvalue};
    x = double(x); y = double(y);
    % Projection conversion
    [x,y] = m_ll2xy(x,y);   
    N_now = N{zvalue};
    NNan = find(~isnan(x));

    if length(NNan) < Min_length; continue; end
    
    % Make the interpolant using scatteredInterpolant natural
    F = scatteredInterpolant(x(NNan),y(NNan),N_now(NNan),'natural','nearest');
    
    % Get data where depth is equal to rounded zvalue
    %I = find( contour_int*int32(B/contour_int) == zcontour(zvalue));
    % Get data where depth is larger than current contour
%     I = find( B >= zcontour(zvalue) & B < zcontour(zvalue+1));
%     Nb = F(xx(I),yy(I));
%     J = find( );
%     Nb_lw = F(xx(J),yy(J));
%     % Interpolate data to get Nb 
%     Nb(I) = F(xx(I),yy(I));
    J = find( B > zcontour(max(1,zvalue-1)));
    % Get data where depth is larger than or equal to round zvalue
    %J = find( contour_int*int32(B/contour_int) >= zcontour(zvalue));
    N_interp{zvalue} = NaN(length(B),1);
    N_interp{zvalue}(J) = F(xx(J),yy(J));
    
    % Add data at each contour for averaging later
    %Nlast(J) = F(xx(J),yy(J));
    %Nm(J) = Nm(J) + Nlast(J);
    %Count(J) = Count(J) + 1;
end
%
for zvalue = 1:length(zcontour)-1    
    DZ = zcontour(zvalue+1) - zcontour(zvalue);
    % Test whether we have data above us
    if ~isempty(N_interp{zvalue+1})
        J = find( B > zcontour(zvalue) & B <= zcontour(zvalue+1));
        if ~isempty(J)
            % For Nb do linear interp
            dz = min(DZ,(B(J) - zcontour(zvalue)));
            Nb(J) = N_interp{zvalue}(J).*(DZ-dz)/DZ + ...
                    N_interp{zvalue+1}(J).*dz/DZ;

            % For Nm do the integral sum
            % For depths within the current range
            Weight = 0.5*(zcontour(zvalue)+B(J)-zcontour(1))./B(J);
            Nm(J) = Nm(J) + 0.5*Weight.*(N_interp{zvalue}(J)+Nb(J)).*dz;
        end
        % For depths larger than current range
        J = find( B > zcontour(zvalue+1));
        if ~isempty(J)
            Weight = 0.5*(zcontour(zvalue)+zcontour(zvalue+1)-zcontour(1))./B(J);
            Nm(J) = Nm(J) + 0.5*Weight.*...
                            (N_interp{zvalue}(J)+N_interp{zvalue+1}(J))*DZ; 
        end
    else
        J = find( B > zcontour(zvalue));
        if ~isempty(J)
            % For Nb just set equal to the last available contour value
            dz = B(J) - zcontour(zvalue);
            Nb(J) = N_interp{zvalue}(J);

            % For Nm do the integral sum
            % For depths within the current range
            Weight = 0.5*(zcontour(zvalue)+B(J)-zcontour(1))./B(J);
            Nm(J) = Nm(J) + Weight.*Nb(J).*dz;
        end
        break;
    end
end
% Divide Nm by depth
Nm = Nm./(B-zcontour(1));

% Divide by the count of contours
%Nm(Count > 0) = Nm(Count > 0)./Count(Count > 0);
% Constant extrapolation to depths outside of last contour 
%Nb(Nb == 0 & Count > 0) = Nlast(Nb == 0 & Count > 0);
%EOF
end