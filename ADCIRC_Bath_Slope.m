function [Hx,Hy] = ADCIRC_Bath_Slope( EToV,xx,yy,B,SA)
%ADCIRC_Bath_Slope : Gets the bathymetric slopes Hx and Hy at each node on 
%                    an unstructured ADCIRC grid
%
% Inputs: EToV - n x 3 array of triangular elements (node indices)
%         xx   - n x 1 vector of x points (Cartesian)
%         yy   - n x 1 vector of y points (Cartesian)
%         B    - n x 1 vector of bathymetric depths
%         SA   - search bandwidth
%
% Outputs: Hx  - n x 1 vector of bathymetric slope in x direction
%          Hy  - n x 1 vector of bathymetric slope in y direction
%
% BY William Pringle, Oct 7th 2016

%% Get the element areas
A = polyarea(xx(EToV(:,1:3))',yy(EToV(:,1:3))')';

%% Compute the slopes for each element
Hxe = zeros(size(A)); Hye = zeros(size(A));
for n = 1:length(A)
    % Get x differences
    a = [ xx(EToV(n,3)) - xx(EToV(n,2))
          xx(EToV(n,1)) - xx(EToV(n,3))   
          xx(EToV(n,2)) - xx(EToV(n,1)) ] ;
    
    % Get y differences
    b = [ yy(EToV(n,2)) - yy(EToV(n,3))
          yy(EToV(n,3)) - yy(EToV(n,1))   
          yy(EToV(n,1)) - yy(EToV(n,2)) ];     
    
    % Compute Hx
    Hxe(n) = 0.5 * B(EToV(n,:))'*b;
    
    % Compute Hy 
    Hye(n) = 0.5 * B(EToV(n,:))'*a;
end

Hx = zeros(size(xx)); Hy = zeros(size(xx));
ne = length(A);
%% Get slopes for each node
for i = 1:length(xx)
     if i == 1
        [I,~] = find(EToV(:,1:3) == i);
     else
        % search radius is SA
        ns = max(1,min(I)-SA);
        nee = min(ne,max(I)+SA);
        [I,~] = find(EToV(ns:nee,1:3) == i); 
        I = I + double(ns) - 1;
     end
     if isempty(I)
        % Search failed
        [I,~] = find(EToV(:,1:3) == i); 
     end
     Hx(i) = sum(Hxe(I))/sum(A(I));
     Hy(i) = sum(Hye(I))/sum(A(I));
end

end

