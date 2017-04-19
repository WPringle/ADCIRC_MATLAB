function d = dpoly_fp(p,mdl,pv,F,ub)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% p are the mesh points
% mdl kdtree of segment points
% pv is the polygon of small islands to ignore mesh inside of them
% F is the interpolant for the DEM
% ub is the upper bound of the bathy (e.g. 5 m)
% d is the distance from point p to closest point on the contour ps (d is
% negative if the point gives bathy below the upper bounds and outside the
% polygon of the oceanside mesh)

np = size(p,1) ;
if isempty(gcp('nocreate'))
    % SERIAL
    
    % Get distances to segment
    [~,d] = knnsearch(mdl,p);
    
    % Get the out polygon
    in = ~InPolygon(p(:,1),p(:,2),pv(:,1),pv(:,2));
else
    % PARALLEL
    
    % Get distances to segment
    d = zeros(np,1);
    Pool = gcp(); num_p = Pool.NumWorkers;
    for idx = 1:num_p
        ns = int64((idx-1)*np/num_p)+1;
        ne = int64(idx*np/num_p);
        f(idx) = parfeval(Pool,@knnsearch,2,...
            mdl,p(ns:ne,:));
    end
    for idx = 1:num_p
        [idx_t, ~, d_t] = fetchNext(f); % Get results into a cell array
        ns = int64((idx_t-1)*np/num_p)+1;
        ne = int64(idx_t*np/num_p);
        d(ns:ne) = d_t;
    end
    
    % Get the out polygon
    in = zeros(np,1); 
    for idx = 1:num_p
        ns = int64((idx-1)*np/num_p)+1;
        ne = int64(idx*np/num_p);
        f(idx) = parfeval(Pool,@InPolygon,1,...
            p(ns:ne,1),p(ns:ne,2),pv(:,1),pv(:,2));
    end
    for idx = 1:num_p
        [idx_t, in_t] = fetchNext(f); % Get results into a cell array
        ns = int64((idx_t-1)*np/num_p)+1;
        ne = int64(idx_t*np/num_p);
        in(ns:ne) = ~in_t;
    end
end

% Get interpolant
bathy = F(p);
% make in == 0 if bathy is above upper bounds
in(bathy > ub) = 0;
% d is negative if outside polygon and below upper bounds
d=(-1).^(in).*d;
