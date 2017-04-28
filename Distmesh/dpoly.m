function d =dpoly(p,pv,mdl)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% p are the mesh points
% pv is the bounded polygon
% mdl kdtree of segment points

% d is the distance from point, p to closest point on polygon, ps (d is
% negative if inside the bounded polygon, pv and positive if outside)
% iloc is the indice in polygon for closest point from p

np = size(p,1) ;

if isempty(gcp('nocreate'))
    % SERIAL
    if ~isempty(mdl) 
        [~,d] = knnsearch(mdl,p);
    else
        % make distance very large as must be very far from any boundaries
        d = 1d8*ones(length(p),1);
    end
    in = InPolygon(p(:,1),p(:,2),pv(:,1),pv(:,2));
else
    % PARALLEL
    Pool = gcp(); num_p = Pool.NumWorkers;

    if ~isempty(mdl) 
        % Get distances to segment
        d = zeros(np,1);
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
    else
        % make distance very large as must be very far from any boundaries
        d = 1d8*ones(length(p),1);
    end

    % Get the inpolygon
    in = zeros(length(p),1);
    for idx = 1:num_p
        ns = int64((idx-1)*np/num_p)+1;
        ne = int64(idx*np/num_p);
        f(idx) = parfeval(Pool,@inpoly,1,p(ns:ne,:)',pv');  %InPolygon,1,...
                          %p(ns:ne,1),p(ns:ne,2),pv(:,1),pv(:,2));
    end
    for idx = 1:num_p
        [idx_t, in_t] = fetchNext(f); % Get results into a cell array
        ns = int64((idx_t-1)*np/num_p)+1;
        ne = int64(idx_t*np/num_p);
        in(ns:ne) = in_t;
    end
end
% d is negative if inside polygon
d=(-1).^(in).*d;
