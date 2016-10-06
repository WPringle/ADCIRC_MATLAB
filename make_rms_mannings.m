% Put reefs, mangroves and sediment mannings onto mesh
clearvars; clc; close all
%out filename
filename = 'fort.13.rms_mannings';
N_default  = 0.020;
N_reef     = 0.22;
N_mangrove = 0.40;
% load the grid
load ../IDIOMS_v5.16.mat

% % load the reefs and mangroves data
% mangrove = 'E:\Indian_Ocean\DataPack-14_001_WCMC010_MangrovesUSGS2011_v1_3';
% mangrove = [mangrove '/01_Data/14_001_WCMC010_MangroveUSGS2011_v1_3.shp'];
% reef = 'E:\Indian_Ocean\DataPack-14_001_WCMC008_CoralReef2010_v1_3';
% reef = [reef '/01_Data/14_001_WCMC008_CoralReef2010_v1_3.shp'];
% 
% S_man  = shaperead(mangrove);
% man_bx = zeros(length(S_man),1); man_by = zeros(length(S_man),1);
% for k = 1:length(S_man)
%     man_bx(k) = S_man(k).BoundingBox(1,1);
%     man_by(k) = S_man(k).BoundingBox(1,2);
% end 
% 
% S_reef = shaperead(reef);
% reef_bx = zeros(length(S_reef),1); reef_by = zeros(length(S_reef),1);
% for k = 1:length(S_reef)
%     reef_bx(k) = S_reef(k).BoundingBox(1,1);
%     reef_by(k) = S_reef(k).BoundingBox(1,2);
% end 
% CVH = convhull(VX(:,1),VX(:,2)); % Get the convex hull
% I = find(inpolygon(man_bx,man_by,VX(CVH,1),VX(CVH,2)) == 0);  
% S_man(I) = [];
% I = find(inpolygon(reef_bx,reef_by,VX(CVH,1),VX(CVH,2)) == 0); 
% S_reef(I) = [];
% save('E:\Indian_Ocean\IDIOMS_v5\Reefs_and_Mangroves_loc.mat','S_man','S_reef');
load E:\Indian_Ocean\IDIOMS_v5\Reefs_and_Mangroves_loc.mat

%% analyse reefs
reef_bx = zeros(length(S_reef),1); reef_by = zeros(length(S_reef),1);
for k = 1:length(S_reef)
    % Get centres of the polygons
    reef_bx(k) = mean(S_reef(k).BoundingBox(1:2,1));
    reef_by(k) = mean(S_reef(k).BoundingBox(1:2,2));
end 
IDX = knnsearch(VX,[reef_bx reef_by],'k',100); 
mannings_reef = zeros(length(IDX),1); ns = 1; ne = 0;
for k = 1:length(IDX)
    % Find the current polygon points
    xnow = S_reef(k).X; xnow(isnan(xnow)) = [];
    ynow = S_reef(k).Y; ynow(isnan(ynow)) = [];
    K = convhull(xnow,ynow); % Get convex hull of polygon
    In = find(inpolygon(VX(IDX(k,:),1),VX(IDX(k,:),2),xnow(K),ynow(K))==1);
    if ~isempty(In)
        if length(In) == 100
            IDX_n = knnsearch(VX,[reef_bx(k) reef_by(k)],'k',1000); 
            In = find(inpolygon(VX(IDX_n,1),VX(IDX_n,2),xnow(K),ynow(K))== 1);
        end
        ne = ne + length(In);
        if length(In) >= 100
            mannings_reef(ns:ne) = IDX_n(In);
        elseif length(In) >= 1000
            disp('shitty'); return
        else
            mannings_reef(ns:ne) = IDX(k,In);
        end
        ns = ne + 1;
    end
end
mannings_reef = mannings_reef(1:ne);
mannings_reef = unique(mannings_reef);

%% analyse mangroves
man_bx = zeros(length(S_man),1); man_by = zeros(length(S_man),1);
for k = 1:length(S_man)
    % Get centres of the polygons
    man_bx(k) = mean(S_man(k).BoundingBox(1:2,1));
    man_by(k) = mean(S_man(k).BoundingBox(1:2,2));
end 
IDX = knnsearch(VX,[man_bx man_by],'k',100); 
mannings_man = zeros(length(IDX),1); ns = 1; ne = 0;
for k = 1:length(IDX)
    % Find the current polygon points
    xnow = S_man(k).X; xnow(isnan(xnow)) = [];
    ynow = S_man(k).Y; ynow(isnan(ynow)) = [];
    K = convhull(xnow,ynow); % Get convex hull of polygon
    In = find(inpolygon(VX(IDX(k,:),1),VX(IDX(k,:),2),xnow(K),ynow(K))==1);
    if ~isempty(In)
        if length(In) == 100
            IDX_n = knnsearch(VX,[man_bx(k) man_by(k)],'k',1000); 
            In = find(inpolygon(VX(IDX_n,1),VX(IDX_n,2),xnow(K),ynow(K))== 1);
        end
        ne = ne + length(In);
        if length(In) >= 100
            mannings_man(ns:ne) = IDX_n(In);
        elseif length(In) >= 1000
            disp('shitty'); return
        else
            mannings_man(ns:ne) = IDX(k,In);
        end
        ns = ne + 1;
    end
end
mannings_man = mannings_man(1:ne);
mannings_man = unique(mannings_man);

%% Put the two mannings together
mannings = zeros(length(mannings_man) + length(mannings_reef),2);
mannings(:,1) = [mannings_man; mannings_reef];
mannings(1:length(mannings_man),2) = N_mangrove; 
mannings(length(mannings_man)+1:end,2) = N_reef;
mannings = sortrows(mannings,1);

%% Get sediment type data
load E:\Global_Data\Sediment_Distribution\Seabed_mannings_v1.mat
I = find(lon >= min(VX(:,1)) & lon <= max(VX(:,1)));
J = find(lat >= min(VX(:,2)) & lat <= max(VX(:,2)));
[lx, ly] = ndgrid(lon(I),lat(J));
man = man(J,I)';
% interpolate
F = griddedInterpolant(lx,ly,man);
man_sed = F(VX(:,1),VX(:,2));

% make array with node indices
man_sed = [(1:length(VX))', man_sed];
% eliminate those with values equal to default
min_man = min(man(man > 0));
man_sed = man_sed(man_sed(:,2) ~= N_default & man_sed(:,2) >= min_man,:);

%% Put sediment mannings with mangroves and reefs
mannings = [mannings; man_sed];
% delete sediment values where mangroves and reefs exist
[~,IA,~] = unique(mannings(:,1),'first');
mannings = mannings(IA,:);
% sort by indices
mannings = sortrows(mannings,1);
% Get number of nodes not equal to default
ne = length(mannings);
%% Lets write out the file
fid = fopen(filename,'w');
% Print header
fprintf(fid,'%s \n','Spatial attributes description') ;
fprintf(fid,'%d \n',length(B)) ;   
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%s \n','mannings_n_at_sea_floor') ;
fprintf(fid,'%s \n','time/(length**1/3)') ;
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%f \n',N_default) ;  
fprintf(fid,'%s \n','mannings_n_at_sea_floor') ;
%number of nodes not default
fprintf(fid,'%d \n',ne) ; 
% Print out list of nodes for each
for k = 1:ne
    fprintf(fid,'%d \t %f \n',mannings(k,1),mannings(k,2));
end
fclose(fid);