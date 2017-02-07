function [node, arc_o, polygon] = Read_SMS_Map( finputname )
% Read_SMS_Map: Reads a map created in SMS to ouput the polygon data
% Outputs a structure of the main outer polygon, followed by inner ones

fid = fopen(finputname) ;

header = textscan(fid,'%s%s%s',13) ;
disp(header{1})

% Get NODES
type = textscan(fid,'%s',1);
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
op_num = 0;
for ii = 1:length(arc)
    if strcmp(arc(ii).type,'inner')
        polygon.inner = [polygon.inner; ...
                         node(arc(ii).node(1),:); ...
                         arc(ii).vertices(:,1:2); ...
                         node(arc(ii).node(2),:); ...
                         NaN NaN];
    elseif strcmp(arc(ii).type,'outer')  
        op_num = op_num + 1;
    end
end

% Now make outer polygon by stitching together ocean and mainland
% boundaries in correct order
arc_o  = arc;
polygon.mainland = [];
for op = 1:op_num
    for ii = 1:length(arc)
        if strcmp(arc(ii).type,'outer') 
            if op == 1
                polygon.outer = [node(arc(ii).node(1),:); ...
                                 arc(ii).vertices(:,1:2); ...
                                 node(arc(ii).node(2),:)];
                node_b = node(arc(ii).node(2),:);
                make_mainland
                arc(ii) = [];
                break
            else
                if node(arc(ii).node(1),1) == node_b(1)
                    polygon.outer(end+1:end+arc(ii).v_num+1,:) = [....
                                               arc(ii).vertices(:,1:2); ...
                                               node(arc(ii).node(2),:)];
                    node_b = node(arc(ii).node(2),:);
                    make_mainland
                    arc(ii) = [];
                    break
                elseif node(arc(ii).node(2),1) == node_b(1);
                    polygon.outer(end+1:end+arc(ii).v_num+1,:) = [....
                                       flipud(arc(ii).vertices(:,1:2)); ...
                                               node(arc(ii).node(1),:)];
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
figure(1);
plot(polygon.outer(:,1),polygon.outer(:,2))
hold on
plot(polygon.inner(:,1),polygon.inner(:,2))
plot(polygon.mainland(:,1),polygon.mainland(:,2))
%EOF
end

