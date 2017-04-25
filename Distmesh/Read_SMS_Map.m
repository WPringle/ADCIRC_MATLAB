function [node, arc_o, polygon] = Read_SMS_Map( finputname, plot_on, poly_num )
% Read_SMS_Map: Reads a map created in SMS to ouput the polygon data
% Outputs a structure of the main outer polygon, followed by inner ones

% plot_on :  plot the polygon or not (=1 or =0)
% splice  : number of splice*splice splices to make in grid
% poly_num : select the sub-polygon number seperated by splitter boundaries
%            = 0 for all, > 0 for a selection

fid = fopen(finputname) ;
header = textscan(fid,'%s%s%s\n',1) ;
if poly_num == 0
    % Just read basically the whole map when only one exists
    while ~strcmp(header{1},'NODE')
        header = textscan(fid,'%s%s%s',1) ;
        disp([header{1} header{2} header{3}])
    end
else
    nn = 0;
    while nn < poly_num
        header = textscan(fid,'%s%s%s',1) ;
        if strcmp(header{1},'BEGCOV')
            nn = nn + 1;
        end
    end
    while ~strcmp(header{1},'NODE')
        header = textscan(fid,'%s%s%s',1) ;
        disp([header{1} header{2} header{3}])
    end
end
%
% Get NODES
type = header;
if strcmp(type{1},'COVATTS')
    dummy = textscan(fid,'%s%s',1);
    type = textscan(fid,'%s',1);
end
while strcmp(type{1},'NODE')
    % reading info
    pos = textscan(fid,'%s %f %f %f',1);
    id = textscan(fid,'%s %d',1);
    
    % putting in array
    node(id{2},1:2) = [pos{2} pos{3}];
    
    % reading endings
    end_s = textscan(fid,'%s',1);
    type = textscan(fid,'%s',1);
end

% Get ARC
while strcmp(type{1},'ARC')
    % reading info
    id = textscan(fid,'%s %d',1);
    ele = textscan(fid,'%s %f',1);
    nodes = textscan(fid,'%s %d %d',1);
    v_num = textscan(fid,'%s %d',1);
    vertices = textscan(fid,'%f %f %f',v_num{2});
    
    % Put in arc struc
    arc(id{2}).node = [nodes{2} nodes{3}];
    arc(id{2}).v_num = v_num{2};
    arc(id{2}).vertices = [vertices{1} vertices{2} vertices{3}];
    
    % reading miscellaneous
    stuff = textscan(fid,'%s %d',3);
    merg = textscan(fid,'%s %d %d %d %d',1);
    adc = textscan(fid,'%s %d %d %d %d %d',1);
    
    % Put in type
    if adc{3} == 0 
        arc(id{2}).type = 'outer';
        if adc{2} == 2
            arc(id{2}).subtype = 'mainland';
        elseif adc{2} == 1
            arc(id{2}).subtype = 'ocean';
        elseif adc{2} == 0
            arc(id{2}).subtype = 'splitter';
        end
    else
        arc(id{2}).type = 'inner';
        arc(id{2}).subtype = 'island';
    end
    
    % reading endings
    end_s = textscan(fid,'%s',2);
    type = textscan(fid,'%s',1);
end

%% Make polygons
% Make inner polygon and count number of ocean and mainland boundaries
polygon.inner = [];
op_num = 0; sp_num = 0;
for ii = 1:length(arc)
    if strcmp(arc(ii).type,'inner')
        new_island = [node(arc(ii).node(1),:); ...
                         arc(ii).vertices(:,1:2); ...
                         node(arc(ii).node(2),:)];
        cw = ispolycw(new_island(:,1), new_island(:,2));
        if cw == 1
            polygon.inner = [polygon.inner; ...
                             new_island;...
                             NaN NaN];
        else
            polygon.inner = [polygon.inner; ...
                             flipud(new_island);...
                             NaN NaN];    
        end
    elseif strcmp(arc(ii).type,'outer')
        op_num = op_num + 1;
        if  strcmp(arc(ii).subtype,'splitter')
            sp_num = sp_num + 1;
        end
    end
end

% Now make outer polygon by stitching together ocean and mainland
% boundaries in correct order
arc_o  = arc; 
polygon.mainland = []; polygon.outer = [];
for op = 1:op_num
    for ii = 1:length(arc)
        if strcmp(arc(ii).type,'outer')
            if op == 1 
                polygon.outer = [polygon.outer;
                                 node(arc(ii).node(1),:); ...
                                 arc(ii).vertices(:,1:2); ...
                                 node(arc(ii).node(2),:)];
                node_b = node(arc(ii).node(2),:);
                make_mainland
                   arc(ii) = [];
                break
            else
                if node(arc(ii).node(1),1) == node_b(1)
                    polytemp = [polygon.outer;
                                arc(ii).vertices(:,1:2); ...
                                node(arc(ii).node(2),:)];
                    polygon.outer = polytemp;
                    node_b = node(arc(ii).node(2),:);
                    make_mainland
                    arc(ii) = [];
                    break
                elseif node(arc(ii).node(2),1) == node_b(1);
                    polytemp = [polygon.outer;
                                flipud(arc(ii).vertices(:,1:2)); ...
                                node(arc(ii).node(1),:)];
                    polygon.outer = polytemp;
                    node_b = node(arc(ii).node(1),:);
                    make_mainland
                    arc(ii) = [];
                   break
                end
            end
        end
    end
end

    function  make_mainland
       if strcmp(arc(ii).subtype,'mainland')
           polygon.mainland(end+1:end+arc(ii).v_num+3,:) = ...
                            [node(arc(ii).node(1),:); ...
                             arc(ii).vertices(:,1:2); ...
                             node(arc(ii).node(2),:);
                             NaN NaN]; 
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

%     function acw = check_anti_cw 
%         acw_x = 0; acw_y = 0;
%         if polygon.outer(end,1) > polygon.outer(end-1,1) 
%             % If come from left need to go up 
%             if polytemp(end,2) > polygon.outer(end,2)
%                 acw_x = 1; 
%             elseif polytemp(end,2) > polygon.outer(end,2)
%                 acw = 1; return
%             end
%         elseif polygon.outer(end,1) < polygon.outer(end-1,1) 
%             % If come from right need to go down
%             if polytemp(end,2) < polygon.outer(end,2)
%                 acw_x = 1; 
%             elseif polytemp(end,2) == polygon.outer(end,2)
%                 acw = 1; return
%             end
%         end
%         if polygon.outer(end,2) > polygon.outer(end-1,2) 
%             % If come from down need to go left
%             if polytemp(end,1) < polygon.outer(end,1)
%                 acw_y = 1;
%             elseif polytemp(end,1) == polygon.outer(end,1)
%                 acw = 1; return
%             end
%         elseif polygon.outer(end,2) < polygon.outer(end-1,2) 
%             % If come from up need to go right
%             if polytemp(end,1) > polygon.outer(end,1)
%                 acw_y = 1; 
%             elseif polytemp(end,1) == polygon.outer(end,1)
%                 acw = 1; return
%             end
%         end
%         dis_x = polytemp(end,2) - polygon.outer(end,2);
%         dis_y = polytemp(end,1) - polygon.outer(end,1);
%         if acw_x == 1 && acw_y == 1
%             % Both is anti-clockwise so must be anti-clockwise
%             acw = 1;
%         elseif acw_x == 0 && acw_y == 0
%             % Both is clockwise so must be clockwise
%             acw = 0;
%         else 
%             % Have some discrepancy
%             if abs(dis_x) > abs(dis_y)
%                 acw = acw_x;
%             elseif abs(dis_y) > abs(dis_x)
%                 acw = acw_y;
%             end
%         end
%     end
% 
% if splice > 0
%     minp = min(polygon.outer); maxp = max(polygon.outer); 
%     dif_x = maxp(1) - minp(1);
%     dif_y = maxp(2) - minp(2);
%     n = 0;
%     for s_y = splice:-1:1
%         if s_y < splice
%             edge(4) = minp(2) + s_y*dif_y/splice ;
%         else
%             edge(4) = maxp(2);
%         end
%         edge(2) = minp(2) + (s_y-1)*dif_y/splice;
%         for s_x = 1:splice
%             if s_x < splice
%                 edge(3) = minp(1) + s_x*dif_x/splice ;
%             else
%                 edge(3) = maxp(1);
%             end
%             edge(1) = minp(1) + (s_x-1)*dif_x/splice;
%             n = n + 1;
%             if n ~= poly_num; continue; end
%             I = find(polygon.outer(:,1) >= edge(1) & ...
%                      polygon.outer(:,1) <= edge(3) & ...
%                      polygon.outer(:,2) >= edge(2) & ...
%                      polygon.outer(:,2) <= edge(4));  
%             if isempty(I); 
%                 polygon = []; 
%             else
%                 end_edge = find_edge(polygon.outer(I(end),:));
%                 start_edge = find_edge(polygon.outer(I(1),:));
%                 ptemp = [polygon.outer(I,:); polygon.outer(I(1),:)];
%                 if end_edge == start_edge && ~ispolycw(ptemp(:,1),ptemp(:,2))
%                     polygon.outer = ptemp;
%                 else
%                     if start_edge > end_edge 
%                         ed = end_edge:start_edge; 
%                     else
%                         ed = [end_edge:4 1:start_edge]; 
%                     end
%                     polygon.outer = polygon.outer(I,:);
%                     ccc = 0;
%                     for dd = ed
%                         ccc = ccc + 1;
%                         up = dd + 1; if up > 4; up = 1; end
%                         if dd == ed(end) && ccc > 1
%                             if dd == 2
%                                polygon.outer = [polygon.outer;  ...
%                                     min(polygon.outer(1,1),edge(up)),edge(dd) ];
%                             elseif dd == 4
%                                polygon.outer = [polygon.outer; ...
%                                     max(polygon.outer(1,1),edge(up)),edge(dd) ];
%                             elseif dd == 3    
%                                polygon.outer = [polygon.outer; edge(dd) ...
%                                         min(polygon.outer(1,2),edge(up)) ];
%                             elseif dd == 1
%                                 polygon.outer = [polygon.outer; edge(dd) ...
%                                         max(polygon.outer(1,2),edge(up)) ];
%                             end
%                         else
%                             if mod(dd,2) == 0
%                                 polygon.outer = [polygon.outer;  ...
%                                                  edge(up),edge(dd) ]; 
%                             else
%                                 polygon.outer = [polygon.outer;  ...
%                                                  edge(dd),edge(up) ]; 
%                             end
%                         end
%                     end    
%                 end
%                 % Get mainland and islands inside
%                 in = InPolygon(polygon.inner(:,1),polygon.inner(:,2),...
%                                        polygon.outer(:,1),polygon.outer(:,2));
%                 cn = 0;
%                 for ii = 1:length(in)
%                     if in(ii) == 1
%                         cn = cn + 1;
%                         inner_n(cn,:) = polygon.inner(ii,:);
%                     elseif ii > 1
%                         if in(ii) == 0 && in(ii-1) == 1
%                            cn = cn + 1;
%                            inner_n(cn,:) = [NaN, NaN];
%                         end
%                     end
%                 end
%                 if exist('inner_n','var')
%                     polygon.inner = inner_n;
%                 else
%                     polygon.inner = [];
%                 end
%             end
%             break;
%         end         
%     end    
% end
% 
%     function edge1 = find_edge(x)
%         if abs(x(1) - edge(1)) < dif_x*0.01
%             edge1 = 1;
%         elseif abs(x(1) - edge(3)) < dif_x*0.01
%             edge1 = 3;     
%         elseif abs(x(2) - edge(2)) < dif_y*0.01
%             edge1 = 2;
%         elseif abs(x(2) - edge(4)) < dif_y*0.01
%             edge1 = 4;
%         end
%     end

%EOF
end

