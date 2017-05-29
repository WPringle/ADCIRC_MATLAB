function polygon = Read_shapefile( finputname, bbox, min_length, ...
                                   h0, plot_on, polygon )
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
    % Get rid of unwanted components
    F = fieldnames(S);
    D = struct2cell(S);
    S = cell2struct(D(3:4,:), F(3:4));
    
    % Get only the polygons which fit within bbox & and are larger than
    % specified number of points
    nn = 0; I = []; L = 0;
    for s = 1:length(S)
        x_n = S(s).X; y_n = S(s).Y; 
        x_n = x_n(~isnan(x_n)); y_n = y_n(~isnan(y_n));
        m_d = mean(abs(diff([x_n; y_n],[],2)),2); m_d = norm(m_d,2);
        % Ignore small length shapes
        if length(x_n) < h0/m_d*min_length; continue; end
        if any(x_n > bbox(1,1)) && any(x_n < bbox(1,2)) && ...
           any(y_n > bbox(2,1)) && any(y_n < bbox(2,2))   
           % Make sure we also ignore the shapes where the length of 
           % the array thats within bbox is small
           xtemp = get_x_y_in_bbox(S(s),bbox);
           if length(xtemp) < h0/m_d*min_length; continue; end
           nn = nn + 1;
           I(nn) = s;
           L = length(xtemp) + L;
        end
    end
    if ~isempty(S)
        if exist('SG','var')
            % Keep the following polygons
            SG = [SG; S(I)];
        else
            SG = S(I);
        end
    end
end

if ~isempty(polygon)
    open = 0;
else
    polygon.outer = [bbox(1,1) bbox(2,1);
                     bbox(1,1) bbox(2,2);
                     bbox(1,2) bbox(2,2);
                     bbox(1,2) bbox(2,1);
                     bbox(1,1) bbox(2,1)];
    open = 1;
end

%% Find the polygons which are wholey inside the bbox.. set as islands
polygon.inner = NaN(L+length(SG)-1,2); 
polygon.mainland = [];
lb_num = 0; ns = 1; ne = 0;
for s = 1:length(SG)
    x_n = SG(s).X; y_n = SG(s).Y;
    %x_n = x_n(~isnan(x_n)); y_n = y_n(~isnan(y_n));
    if min_length > 0
        if exist('inpoly','file') == 2
            % m-file version
            In = inpoly([x_n(~isnan(x_n))',y_n(~isnan(y_n))'],polygon.outer);
        elseif exist('inpoly','file') == 3
            % Mex version
            In = inpoly([x_n(~isnan(x_n));y_n(~isnan(y_n))],polygon.outer');
        end
    else
        In = ones(length(x_n(~isnan(x_n))),1);
    end
    if length(find(In == 1)) == length(x_n(~isnan(x_n)))
    %if min(x_n) > bbox(1,1) && max(x_n) < bbox(1,2) && ...
    %   min(y_n) > bbox(2,1) && max(y_n) < bbox(2,2) 
        new_island = [x_n' y_n'];
        if min_length > 0
            cw = ispolycw(new_island(:,1), new_island(:,2));
        else
            cw = 1;
        end
        ne = ne + size(new_island,1);
        if ne > size(polygon.inner,1)
           polygon.inner = [polygon.inner; NaN(length(SG),2)];
        end
        if cw == 1
            polygon.inner(ns:ne,:) = new_island;
        else
            polygon.inner(ns:ne,:) = flipud(new_island);    
        end
        ns = ne + 1;
    else
        if open == 0
            x_n = x_n(~isnan(x_n)); y_n = y_n(~isnan(y_n));
            new_main = [x_n(In == 1)' y_n(In == 1)'];
            if ~isempty(new_main)
                polygon.mainland = [polygon.mainland; new_main; NaN NaN]; 
            end
        else
            lb_num = lb_num + 1;
            lb_v(lb_num) = s;
        end
    end
end
%polygon.inner(polygon.inner(:,1) == 0 & polygon.inner(:,2) == 0,:) = [];
% Now make outer polygon by stitching together ocean and mainland
% boundaries in correct order
% if open = 0 just make mainland
if open == 0
    if plot_on == 1 && ~isempty(polygon)
        figure(1);
        hold on
        if ~isempty(polygon.inner)
            plot(polygon.inner(:,1),polygon.inner(:,2))
        end
        if ~isempty(polygon.mainland)
            plot(polygon.mainland(:,1),polygon.mainland(:,2))
        end
    end
    return;
end
polygon.outer = []; 
end_e = 2;
% Get the starting lb for the longest shape within bbox
mL = 0;
for lbt = 1:length(lb_v)
    X_n = get_x_y_in_bbox(SG(lb_v(lbt)),bbox);
    if length(X_n) > mL
        lb = lbt;
    end
end
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
                segment = add_corner(polygon.outer(end,:),e,bbox,h0);
                polygon.outer(end+1:end+length(segment),:) = segment;
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
    if length(I) < 2; return; end
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

function corner = add_corner(poly_end,e,bbox,h0)
% find the current edge number, = 1 top, = 2 right, =3 bottom, = 4 left
    if e == 1
        num = abs(poly_end(1) - bbox(1,2))/h0;
        corner(:,1) = linspace(poly_end(1),bbox(1,2),num);
        %corner(:,1) = poly_end(1)-h0:-h0:bbox(1,2); 
        corner(end,1) = bbox(1,2); corner(:,2) = bbox(2,2);
        %corner = [bbox(1,2) bbox(2,2)];
    elseif e == 2
        %corner(:,2) = poly_end(2)+h0:h0:bbox(2,1); 
        num = abs(poly_end(2) - bbox(2,1))/h0;
        corner(:,2) = linspace(poly_end(2),bbox(2,1),num); 
        corner(end,2) = bbox(2,1); corner(:,1) = bbox(1,2); 
        %corner = [bbox(1,2) bbox(2,1)];
    elseif e == 3
        num = abs(poly_end(1) - bbox(1,1))/h0;
        corner(:,1) = linspace(poly_end(1),bbox(1,1),num); 
        %corner(:,1) = poly_end(1)+h0:h0:bbox(1,1); 
        corner(end,1) = bbox(1,1);
        corner(:,2) = bbox(2,1);
        %corner = [bbox(1,1) bbox(2,1)];
    elseif e == 4
        num = abs(poly_end(2) - bbox(2,2))/h0;
        corner(:,2) = linspace(poly_end(2),bbox(2,2),num); 
        %corner(:,2) = poly_end(2)-h0:-h0:bbox(2,2); 
        corner(end,2) = bbox(2,2); corner(:,1) = bbox(1,1); 
        %corner = [bbox(1,1) bbox(2,2)]; 
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

