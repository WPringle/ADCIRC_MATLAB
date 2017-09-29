function polygon_struct = Read_shapefile( finputname, polygon, bbox, ...
    h0, plot_on )
% Read_shapefile: Reads a shapefile or a NaN-delimited vector
% containing polygons and/or segments in the the desired region 
% of interest. Classifies the vector data as either a
% mainland, island, or outer boundary. The program will automatically
% trim islands that are have an area smaller than 4*h0^2 and will in gaps
% in the vector that are larger than h0/2. This is necessary for the
% use of boundary rejection method in dpoly. 

% INPUTS: 
% finputname : file name(s) of the shapefile listed in a cell array
% polygon    : a NaN-delimited vector of polygons and/or segments.
% bbox       : the bounding box that we want to extract out
% h0         : minimum edge length in region of interest. 
% plot_on    : plot the final polygon or not (=1 or =0)
%
% OUTPUTS: 
% polygon_struct    : a structure containing the vector features identified as
%              either islands, mainland, or outer.
% Written by William Pringle and Keith Roberts, CHL,UND, 2017
%% Loop over all the filenames and get the shapefile within bbox
SG = [];
if(size(finputname,1)~=0)
    for fname = finputname
        % Read the structure
        S = shaperead(fname{1},'BoundingBox',bbox');
        
        % Get rid of unwanted components
        F = fieldnames(S);
        D = struct2cell(S);
        S = cell2struct(D(3:4,:), F(3:4));
        
        if ~isempty(S)
            % Keep the following polygons
            SG = [SG; S];
        end
    end
else
    count = 1;
    j=1;
    for i = 1 : length(polygon)
        if(isnan(polygon(i,1))==1)
            count = count + 1; j=1;
            continue 
        end
        SG(count).X(j) = polygon(i,1);
        SG(count).Y(j) = polygon(i,2);
        j=j+1; 
    end
end
% If we don't have an outer polygon already then make it by bbox
polygon_struct.outer = [bbox(1,1) bbox(2,1);
    bbox(1,1) bbox(2,2);
    bbox(1,2) bbox(2,2);
    bbox(1,2) bbox(2,1);
    bbox(1,1) bbox(2,1)];
% Densify the outer polygon (fills gaps larger than half min edgelength).
[latout,lonout] = interpm(polygon_struct.outer(:,2),...
    polygon_struct.outer(:,1),h0/2);

polygon_struct.outer = [];
polygon_struct.outer(:,1) = lonout;
polygon_struct.outer(:,2) = latout;

    %% Find whether the polygons are wholey inside the bbox..
    %% Set as islands or mainland
    polygon_struct.inner = [];
    polygon_struct.mainland = [];
    edges = Get_poly_edges( [polygon_struct.outer;NaN NaN] );

    for s = 1:length(SG)
        % Get current polygon
        x_n = SG(s).X; y_n = SG(s).Y;
        x_n = x_n(~isnan(x_n)); y_n = y_n(~isnan(y_n));
        % Check proportion of polygon that is within bbox
        In = inpoly([x_n',y_n'],[polygon_struct.outer;NaN NaN],edges);
        
        % lets calculate the area of the
        % feature using the shoelace algorithm and decided whether to keep or
        % not based on area.
        if(x_n(end)==x_n(1))
            area = shoelace(x_n',y_n');
        else
            area = 999; % not a polygon
        end
        if(area < 4*h0^2) % too small, then don't consider it.
            continue;
        elseif length(find(In == 1)) == length(x_n)
            % Wholey inside box, set as island
            new_island = [SG(s).X' SG(s).Y'];
            polygon_struct.inner = [polygon_struct.inner; new_island];
        else
            % Partially inside box, set as mainland
            new_main = [SG(s).X' SG(s).Y'];
            polygon_struct.mainland = [polygon_struct.mainland; new_main];
        end
    end
    % Add mainland to outer polygon to get full outer polygon
    polygon_struct.outer = [polygon_struct.outer; NaN NaN; polygon_struct.mainland];
    
    %% Plot the map
    if plot_on >= 1 && ~isempty(polygon_struct)
        figure(1);
        hold on
        plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2))
        if ~isempty(polygon_struct.inner)
            plot(polygon_struct.inner(:,1),polygon_struct.inner(:,2))
        end
        if ~isempty(polygon_struct.mainland)
            plot(polygon_struct.mainland(:,1),polygon_struct.mainland(:,2))
        end
    end
    %EOF
end
