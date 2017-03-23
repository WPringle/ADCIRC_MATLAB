function polygon = Read_shapefile( finputname, bbox, min_length, plot_on )
% Read_shapefile: Reads a shapefile polygon on the coastline, extracting
% the desired area out, and making the open ocean boundaries automatically

% finputname : file name(s) of the shapefile
% bbox    : the bounding box that we want to extract out
% min_length : minimum length of points in polygon to accept 
% plot_on :  plot the final polygon or not (=1 or =0)

%% Loop over all the filenames
for fname = finputname
    % Read the structure
    S = shaperead(fname{1});
    % Get only the polygons which fit within bbox & and are larger than
    % specified number of points
    nn = 0; I = [];
    for s = 1:length(S)
        x_n = S(s).X; y_n = S(s).Y;
        if length(x_n) < min_length; continue; end
        if any(x_n > bbox(1,1)) && any(x_n < bbox(1,2)) && ...
           any(y_n > bbox(2,1)) && any(y_n < bbox(2,2))   
           nn = nn + 1;
           I(nn) = s;
        end
    end
    if exist('SG','var')
        % Keep the following polygons
        SG = [SG; S(I)];
    else
        SG = S(I);
    end
end

%% Find the polygons which are wholey inside the bbox.. set as islands
polygon.inner = [];
lb_num = 0; 
for s = 1:length(SG)
    x_n = SG(s).X; y_n = SG(s).Y;
    if min(x_n) > bbox(1,1) && max(x_n) < bbox(1,2) && ...
       min(y_n) > bbox(2,1) && max(y_n) < bbox(2,2) 
        new_island = [x_n' y_n'];
        cw = ispolycw(new_island(:,1), new_island(:,2));
        if cw == 1
            polygon.inner = [polygon.inner; ...
                             new_island];
        else
            polygon.inner = [polygon.inner; ...
                             flipud(new_island)];    
        end
    else
        lb_num = lb_num + 1;
        lb_v(lb_num) = s;             
    end
end

% Now make outer polygon by stitching together ocean and mainland
% boundaries in correct order
polygon.mainland = []; polygon.outer = []; lb = 1; end_e = 2;
while ~isempty(lb_v)
    X_n = get_x_y_in_bbox(SG(lb_v(lb)),bbox);
    if isempty(X_n)
        lb_v(lb) = [];
        continue; 
    end
    if end_e == 1
        % add to the mainland boundary
        polygon.mainland = [polygon.mainland; X_n; NaN NaN]; 
        % add to the outer polygon                
        polygon.outer = [polygon.outer; X_n];   
    elseif end_e == 2
        % add to the mainland boundary
        polygon.mainland = [polygon.mainland; flipud(X_n); NaN NaN]; 
        % add to the outer polygon                
        polygon.outer = [polygon.outer; flipud(X_n)];
    end
    % delete this boundary
    lb_v(lb) = [];             
    % now we need to draw the open boundary
    % find the current edge number, = 1 top, = 2 right, =3 bottom, = 4 left
    edge = find_edge_num(polygon.outer(end,:),bbox);
    e_vec = edge:4;
    if edge <= 4; e_vec = [e_vec, 1:edge-1]; end
    % now find closest node near edge going around in clockwise direction
    dist = 9e6;
    for e = e_vec
        found = 0; delete = 0;
        for llb = 1:length(lb_v)
            X_n = get_x_y_in_bbox(SG(lb_v(llb)),bbox);
            if isempty(X_n)
                delete = llb;
                continue; 
            end
            se = find_edge_num(X_n(1,:),bbox);
            ee = find_edge_num(X_n(end,:),bbox);
            if se == e
               found = 1;
               d_n = distance(X_n(1,2),X_n(1,1),...
                              polygon.outer(end,2),polygon.outer(end,1));
               if d_n < dist
                    lb = llb; end_e = 1; dist = d_n;
               end
            end
            if ee == e
               found = 1;
               d_n = distance(X_n(end,2),X_n(end,1),...
                              polygon.outer(end,2),polygon.outer(end,1));
               if d_n < dist
                    lb = llb; end_e = 2; dist = d_n;
               end
            end
        end
        if delete > 0
            lb_v(delete) = [];
        end
        if found
            break; 
        else
            ee = find_edge_num(polygon.outer(1,:),bbox);
            if isempty(lb_v) && ee == e
                % set the first point to close polygon
                polygon.outer(end+1,:) = polygon.outer(1,:);
                break;
            else
                % start a new edge so need to add the the next corner in
                polygon.outer(end+1,:) = add_corner(e,bbox);
            end
        end
    end
end  
if isempty(polygon.outer)
   % Just join the four points
   polygon.outer = [bbox(1,2) bbox(2,2);
              bbox(1,2) bbox(2,1);
              bbox(1,1) bbox(2,1);
              bbox(1,1) bbox(2,2);
              bbox(1,2) bbox(2,2)]; 
end

%% Plot the map
if plot_on == 1 && ~isempty(polygon)
    figure(1);
    hold on
    plot(polygon.outer(:,1),polygon.outer(:,2))
    if ~isempty(polygon.inner)
        plot(polygon.inner(:,1),polygon.inner(:,2))
    end
    if ~isempty(polygon.mainland)
        plot(polygon.mainland(:,1),polygon.mainland(:,2))
    end
end
%EOF
end

function X_n = get_x_y_in_bbox(SG,bbox)
    X_n(:,1) = SG.X'; X_n(:,2) = SG.Y';
    I = find(X_n(:,1) > bbox(1,1) & X_n(:,1) < bbox(1,2) & ...
             X_n(:,2) > bbox(2,1) & X_n(:,2) < bbox(2,2));
    X_n = X_n(I,:);
    if isempty(I); return; end
    % Find the points that are close to an edge and re-sort if necessary
    is = [];
    m_d = mean(abs(diff(X_n))); m_d = norm(m_d,2);
    for i = 1:length(X_n)
        [~,d] = find_edge_num(X_n(i,:),bbox);
        if d < 2*m_d
           is = [is i];
        end
    end
    if ~isempty(is)
        % only re-sort if all points are not equal to 1 and X_n
        if ~any(is == 1) && ~any(is == length(X_n))
            X_n = X_n([is:end 1:is-1],:);
        end
    end
end

function [edge,d] = find_edge_num(node,bbox)
% find the current edge number, = 1 top, = 2 right, =3 bottom, = 4 left
	[d,edge] = min([bbox(2,2) - node(2);...
                   bbox(1,2) - node(1); ...
                   node(2) - bbox(2,1);...
                   node(1) - bbox(1,1)]);
end

function corner = add_corner(e,bbox)
% find the current edge number, = 1 top, = 2 right, =3 bottom, = 4 left
    if e == 1
        corner = [bbox(1,2) bbox(2,2)];
    elseif e == 2
        corner = [bbox(1,2) bbox(2,1)];
    elseif e == 3
        corner = [bbox(1,1) bbox(2,1)];
    elseif e == 4
        corner = [bbox(1,1) bbox(2,2)]; 
    end
end

function d = find_distance_to_corner(nodes,bbox)
% set up the corners
    corners = [bbox(1,2) bbox(2,2);
              bbox(1,2) bbox(2,1);
              bbox(1,1) bbox(2,1);
              bbox(1,1) bbox(2,2)]; 
% get nearest distances
    [~,d] = knnsearch(corners,nodes);     
end

