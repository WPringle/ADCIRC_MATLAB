function [sponge,opedat,boudat,pv,B] = get_latlon_for_sponge_zone(...
                                                 f14file,period,frac,write)
% [pvv,idx_g] = get_latlon_for_sponge_zone(f14file,period,frac)
% This function reads the open boundaries of the mesh and gets the nodes in
% the sponge layer as specified by the length, L
% Inputs:
% period = 12.42*3600; % period of the M2 wave
% frac   = 0.1;        % length of the sponge zone/length of M2 wave
% write  = 1 write out spng_lat_lon to file

g = 9.81; % gravity

% Load mesh
if strcmp(f14file(end-2:end),'grd') || strcmp(f14file(end-1:end),'14')
    [~,pv,B,opedat,boudat] = readfort14(f14file);
elseif strcmp(f14file(end-2:end),'mat')
    load(f14file)
    if exist('VX','var')
       pv = VX; clear VX; 
    end
end

% Now loop over boundaries and get the sponge layer points
idx_g = []; 
for op = 1:opedat.nope
    nodes = opedat.nbdv(1:opedat.nvdll(op),op);
    % get mean depth along boundary
    H = mean(B(nodes));
    % Get the length of the sponge (in metres)
    L = frac(op)*period*sqrt(g*H); 
    % Change length to degrees (roughly)
    L = L*1d-5;
    % Find all points L distance from the mean lon or lat of open boundary
    idx = cell(2,1);
    if mean(abs(diff(pv(nodes,1))))/mean(abs(diff(pv(nodes,2)))) < 0.05
       % boundary is essentially parallel to longitude
       sponge(op).orientation = 1;
       for ii = 1:2
            % go backwards and forwards of line and 
            % keep direction of most number of nodes
            xmax = max(pv(nodes,1))*(2-ii) + min(pv(nodes,1))*(ii-1) + (ii-1)*L;
            xmin = max(pv(nodes,1))*(2-ii) + min(pv(nodes,1))*(ii-1) + (ii-2)*L;
            % Get max, min nodes in y
            idx{ii} = find(pv(:,1) <= xmax & pv(:,1) >= xmin & ...
                       pv(:,2) < max(pv(nodes,2)) + 1 & ...
                       pv(:,2) > min(pv(nodes,2)) - 1);   
       end
    elseif mean(abs(diff(pv(nodes,2))))/mean(abs(diff(pv(nodes,1)))) < 0.05
        % boundary is essentially parallel to latitude
        sponge(op).orientation = 2;
        for ii = 1:2
            % go backwards and forwards of line and 
            % keep direction of most number of nodes
            ymax = max(pv(nodes,2))*(2-ii) + min(pv(nodes,2))*(ii-1) + (ii-1)*L;
            ymin = max(pv(nodes,2))*(2-ii) + min(pv(nodes,2))*(ii-1) + (ii-2)*L;
            idx{ii} = find(pv(:,2) <= ymax & pv(:,2) >= ymin & ...
                       pv(:,1) < max(pv(nodes,1)) + 1 & ...
                       pv(:,1) > min(pv(nodes,1)) - 1);   
        end
    else
        % We have a curved boundary
        sponge(op).orientation = 3;
        % Get the nodes in the trivial box from OB
        xmax = max(pv(nodes,1)) + L; xmin = min(pv(nodes,1)) - L;
        ymax = max(pv(nodes,2)) + L; ymin = min(pv(nodes,2)) - L;
        I = find(pv(:,2) <= ymax & pv(:,2) >= ymin & ...
                 pv(:,1) <= xmax & pv(:,1) >= xmin);   
        % Do dsegment on all the I nodes to get distance to the OB
        d = dsegment(pv(I,:),pv(nodes,:));
        d = min(d,[],2);
        idx{1} = (I(d < L));
    end
    [~,maxloc] = max(cellfun('length', idx));
    if maxloc == 1
         sponge(op).orientation = sponge(op).orientation*-1;
    end
    idx = idx{maxloc};
    idx_g = [idx_g; idx]; 
    sponge(op).idx = idx;
    sponge(op).pv  = pv(idx,:);
    sponge(op).H   = H;
    sponge(op).L   = L;
end
if write
    % write lat and lon
    dlmwrite('spng_lat_lon',fliplr(pv(idx_g,:)),'precision',7);
end
%EOF
end