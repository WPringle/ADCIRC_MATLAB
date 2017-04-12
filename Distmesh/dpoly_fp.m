function [d,iloc]= dpoly_fp(p,ps,pv,F,bounds)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% p are the mesh points
% ps is the segment for the contours in floodplain (e.g. 0 m and 5 m
% contours) 
% pv is the polygon of small islands to ignore mesh inside of them
% F is the interpolant for the DEM
% bounds(1) is the lower bound of the bathy (e.g. 0 m)
% bounds(2) is the upper bound of the bathy (e.g. 5 m)

% d is the distance from point p to closest point on the contour ps (d is
% negative if the point gives bathy within the bounds specified)

np=size(p,1) ;
if ~isempty(ps)
    %nvs=size(pv,1)-1;
    nvs=size(ps,1)-1;
    memneed = np*nvs*8/10^9  ; % assume real number  
    if isempty(gcp('nocreate'))
        %disp('dsegment serial')
        maxmem = 4; % 4 Gb
        if ( memneed < maxmem ) % 8 Gb
           %ds=dsegment(p,pv) ;
           ds=dsegment(p,ps) ;
           [d,~]=min(ds,[],2);
        else
           csz = ceil(maxmem*10^9/(8*nvs)) ;

           ires = mod( np, csz) ;
           if ( ires )
               ichunck = [0:csz:np np] ;
           else
               ichunck = [0:csz:np] ; 
           end

           d = zeros(np,1) ;
           iloc = zeros(np,1) ; 
           nch = length(ichunck) ;
           for ic = 1: nch - 1
               % ic
               ibeg = ichunck(ic) + 1 ;
               iend = ichunck(ic+1)   ;

               pp = p(ibeg:iend,:) ; 

               %dsc = dsegment(pp,pv) ;
               dsc = dsegment(pp,ps) ;
               [d(ibeg:iend),iloc(ibeg:iend)] = min( dsc, [], 2) ;  
           end
        end
    else
        %disp('dsegment parallel')
        d = zeros(length(p),1); iloc = zeros(length(p),1);
        parfor ii = 1:np
            ds              = dsegment(p(ii,:),ps);
           [d(ii),iloc(ii)] = min(ds);
        end
    end
else
    % make distance very large as must be very far from any boundaries
    iloc = zeros(length(p),1);
    d = 1d8*ones(length(p),1);
end
% Check if in the bounds or not
% Check to see if in small polygon
in = ~InPolygon(p(:,1),p(:,2),pv(:,1),pv(:,2));
% Get interpolant
bathy = F(p);
% make in == 1 if depth is within the bounds
in(bathy > bounds(1) & bathy < bounds(2)) = 1;
d=(-1).^(in).*d;
%disp('dpoly out')
%
