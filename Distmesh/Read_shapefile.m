function polygon = Read_shapefile( finputname, bbox, min_length, plot_on )
% Read_shapefile: Reads a shapefile polygon on the coastline, extracting
% the desired area out, and making the open ocean boundaries automatically

% finputname : file name(s) of the shapefile
% bbox    : the boundaing box that we want to extract out
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
polygon.mainland = []; polygon.outer = [];
for s = 1:lb_num
    x_n = SG(lb_v(s)).X; y_n = SG(lb_v(s)).Y;
    I = find(x_n > bbox(1,1) & x_n < bbox(1,2) & ...
             y_n > bbox(2,1) & y_n < bbox(2,2));      
    polygon.mainland = [polygon.mainland;
                        x_n(I)' y_n(I)';
                        NaN NaN];               
    polygon.outer = [polygon.outer;
                     x_n(I)' y_n(I)'];
end
% making corner array
corners = [bbox(1,1) bbox(2,1);
           bbox(1,2) bbox(2,1);
           bbox(1,2) bbox(2,2);
           bbox(1,1) bbox(2,2)];
% Get starting & ending corner
c_end  = knnsearch(corners,polygon.outer(end,:));
c_start = knnsearch(corners,polygon.outer(1,:));
c_vec = c_end:4;
if c_end < 4; c_vec = [c_vec, 1:c_end-1]; end
% loop over corners in acw direction making open boundary
for c = 1:4
    if c_vec(c) == 1 || c_vec(c) == 4
        polygon.outer(end+1,:) = [polygon.outer(end,2) corners(c_vec(c),2)];
    elseif c_vec(c) == 2 || c_vec(c) == 3
        polygon.outer(end+1,:) = [corners(c_vec(c),1) polygon.outer(end,2)];
    end
    if c_start == c_vec(c)
      polygon.outer(end+1,:) = polygon.outer(1,:);
      break;
    end
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

