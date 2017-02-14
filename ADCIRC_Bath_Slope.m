function [Hx,Hy] = ADCIRC_Bath_Slope( EToV,xx,yy,B)
%ADCIRC_Bath_Slope : Gets the bathymetric slopes Hx and Hy at each node on 
%                    an unstructured ADCIRC grid
%
% Inputs: EToV - nE x 3 array of triangular elements (node indices)
%         xx   - n x 1 vector of x points (Cartesian)
%         yy   - n x 1 vector of y points (Cartesian)
%         B    - n x 1 vector of bathymetric depths
%
% Outputs: Hx  - n x 1 vector of bathymetric slope in x direction
%          Hy  - n x 1 vector of bathymetric slope in y direction
%
% Created by William Pringle, Oct 7th 2016
% Updated by William Pringle, Dec 2nd 2016 for faster performance
%tic
%% Get the element areas
A = polyarea(xx(EToV(:,1:3))',yy(EToV(:,1:3))')';

%% Compute the slopes for each element and sum over all nodes of that element
Hx = zeros(size(xx)); Hy = zeros(size(xx));
An  = zeros(size(xx));
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
    Hxe = 0.5 * B(EToV(n,:))'*b;
    
    % Compute Hy 
    Hye = 0.5 * B(EToV(n,:))'*a;
    
    % Add in sum for nodes
    Hx(EToV(n,:)) = Hx(EToV(n,:)) + Hxe;
    Hy(EToV(n,:)) = Hy(EToV(n,:)) + Hye;
    An(EToV(n,:)) = An(EToV(n,:)) + A(n);
end
Hx = Hx./An;
Hy = Hy./An;
%toc
end