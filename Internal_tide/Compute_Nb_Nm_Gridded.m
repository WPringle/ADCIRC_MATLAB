function [Nb,Nm,Nmw] = Compute_Nb_Nm_Gridded(lon_M,lat_M,B,zcontour,N,lon_N,lat_N)
                                  %,lon0,lat0)
% Compute_Nb_Nm_Gridded: Compute the buoyancy frequency, N at the seabed 
%                (Nb) and the mean over the depth (Nm) for an unstructured 
%                mesh from grid points of N at specified contours
%
% [Nb,Nm,Nmw] = Compute_Nb_Nm(lon_M,lat_M,B,zcontour,N,lon_N,lat_N,lon0,lat0)
% Input : lon_M    - longitude points of nodes in mesh
%         lat_M    - latitude points of nodes in mesh
%         B        - depths of nodes in mesh
%         zcontour - the contours where we have values of N
%         N        - matrix of N  (lon,lat,z)
%         lon_N    - vector of lon 
%         lat_N    - vector of lat
%         proj     - string that defines the projection to use
%                    e.g. 'Mercator'.., (type m_proj('set') for options)
%
% Output : Nb     - Buoyancy frequency at seabed
%          Nm     - Depth-averaged buoyancy frequency over depth
%          Nmw    - Weighted depth-average buoyancy frequency (linearly
%                   decreasing weights from bottom to surface)
%
% Author: William Pringle, CHL, Notre Dame University
% Created: 2017-9-28

[Lon,Lat] = ndgrid(lon_N,lat_N);
%% Calculation
% initialisation
Nb = zeros(size(B)); Nm = zeros(size(B)); Nmw = zeros(size(B)); 
% do the interpolation onto the mesh
N_interp = cell(length(zcontour),1);
for zvalue = 1:length(zcontour)
    
    % Make the interpolant using griddedInterpolant linear
    F = griddedInterpolant(Lon,Lat,squeeze(N(:,:,zvalue)));
    
    % Find all nodes less than current depth
    J = find( B > zcontour(max(1,zvalue-1)));
    
    % Make the interp cell and interpolate into it
    N_interp{zvalue}    = NaN(size(B));
    N_interp{zvalue}(J) = F(lon_M(J),lat_M(J));
    % Nearest neighbour fill for NaN results
    while ~isempty(find(isnan(N_interp{zvalue}(J)), 1))
        N_interp{zvalue} = fillmissing(N_interp{zvalue},'nearest',1);
        if size(N_interp{zvalue},2) > 1
            N_interp{zvalue} = fillmissing(N_interp{zvalue},'nearest',2);
        end
    end
end
%
for zvalue = 1:length(zcontour)  
    % Test whether we have data above us
    if zvalue < length(zcontour)
        DZ = zcontour(zvalue+1) - zcontour(zvalue);
        J = find( B > zcontour(zvalue) & B <= zcontour(zvalue+1));
        if ~isempty(J)
            % For Nb do linear interp
            dz = min(DZ,(B(J) - zcontour(zvalue)));
            Nb(J) = N_interp{zvalue}(J).*(DZ-dz)/DZ + ...
                    N_interp{zvalue+1}(J).*dz/DZ;

            % For Nm do the integral sum
            % For depths within the current range
            Nm(J) = Nm(J) + 0.5*(N_interp{zvalue}(J)+Nb(J)).*dz;
            Weight = 0.5*(zcontour(zvalue)+B(J))./B(J);
            Nmw(J) = Nmw(J) + 0.5*Weight.*(N_interp{zvalue}(J)+Nb(J)).*dz;
        end
        % For depths larger than current range
        J = find( B > zcontour(zvalue+1));
        if ~isempty(J)
            Nm(J) = Nm(J) + 0.5*(N_interp{zvalue}(J)+N_interp{zvalue+1}(J))*DZ; 
            Weight = 0.5*(zcontour(zvalue+1)+zcontour(zvalue))./B(J);
            Nmw(J) = Nmw(J) + 0.5*Weight.*...
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
            Nm(J) = Nm(J) + Nb(J).*dz;
            Weight = 0.5*(zcontour(zvalue)+B(J))./B(J);
            Nmw(J) = Nmw(J) + Weight.*Nb(J).*dz;
        end
    end
end
% Divide Nm by depth
Nm = Nm./(B-zcontour(1));
Nmw = Nmw./(B-zcontour(1));
Nb(isnan(Nb)) = 0;
Nm(isnan(Nm)) = 0;
Nmw(isnan(Nmw)) = 0;
%EOF
end