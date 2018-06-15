%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read kml file created in google maps to get coordinates of stations 
% (we can move the stations to reasonable locations off land etc in google
% maps before exporting to kml). 
%
% Automatically finds the boundary of your grid and selects only those
% ones that lie within the grid to output to the .csv for copying into
% fort.15 
%
% Requires: kml2struct.m found at https://www.mathworks.com/matlabcentral/
%                                         fileexchange/35642-kml2struct
%           inpoly.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clearvars; clc; close all;

%% Variables to set
kmlfile = {'May14_UpdatedStationPositions.kml'};

% Set the constituents we wanna look at
const = {'M2','S2','N2','K2','K1','O1','P1','Q1'};

% Set output file name
output = 'PRVI_UpdatedStationPositions_NoRepeats.csv';

% Sources in order of priority/reliability
sources = {'Truth_Pelagic','Truth_Shelf','NOAA','JMA','AUS_Tides',...
           'KHAO','GESLA','UHSLC_Fast','Truth_Coast',...
           'Various_Publications','ST727','IHO',};

%% Reading the kml
T = [];
for kml = kmlfile
    kmlS = kml2struct(kml{1});
    % Getting rid of unnecessary shit
    T_n = struct2table(kmlS);
    T_n.Geometry = [];
    T_n.BoundingBox = [];
    %T.Description = [];
    T_n = T_n(:,[1 3 4 2 5]);
    T_n.Source = T_n.StyleUrl;
    T_n.StyleUrl = [];
    T = [T; T_n];
end

%% Convert description to amp and phases
obs = NaN(height(T),length(const)*2);
for i = 1:height(T)                      
    if strcmp(T.Source{i}(end-5:end),'0288D1') % light blue
        T.Source{i} = 'Truth_Pelagic';   
    elseif strcmp(T.Source{i}(end-5:end),'FFEA00') % yellow
        T.Source{i} = 'Truth_Shelf';     
    elseif strcmp(T.Source{i}(end-5:end),'0F9D58') % green
        T.Source{i} = 'Truth_Coast';    
    elseif strcmp(T.Source{i}(end-5:end),'757575') % grey
        T.Source{i} = 'NOAA';  
        % Trim the quotes
        T.Name{i} = T.Name{i}(2:end-1);
        %disp(T.Name{i})
    elseif strcmp(T.Source{i}(end-5:end),'673AB7') % dark purple
        T.Source{i} = 'JMA';            
    elseif strcmp(T.Source{i}(end-5:end),'FF5252') % coral
        T.Source{i} = 'KHAO';           
    elseif strcmp(T.Source{i}(end-5:end),'7CB342') % light green
        T.Source{i} = 'GESLA';           
    elseif strcmp(T.Source{i}(end-5:end),'A52714') % medium red
        T.Source{i} = 'UHSLC_Fast';     
    elseif strcmp(T.Source{i}(end-5:end),'FFD600') % austrlian gold
        T.Source{i} = 'AUS_Tides';     
    elseif strcmp(T.Source{i}(end-5:end),'880E4F') % dark blue
        T.Source{i} = 'ST727';             
    elseif strcmp(T.Source{i}(end-5:end),'3949AB') % dark blue
        T.Source{i} = 'IHO';             
    elseif strcmp(T.Source{i}(end-5:end),'9C27B0') % light purple 
        T.Source{i} = 'Various_Publications'; 
    end
    Descrip = T.Description(i); Descrip = Descrip{1};
    cc = 0;
    for c = const
        cc = cc + 1;
        % Get position in description of the constituent
        Positions = regexp(Descrip,[c{1} '_']);
        Positions = find(~cellfun(@isempty,Positions)); 
        if isempty(Positions); continue; end
        if length(Positions) > 2
            % M2 seasonal shit
            Positions(1:2) = [];
        end
        amp = Descrip(Positions(1));
        s = strfind(amp,' ');
        obs(i,cc*2-1) = str2double(amp{1}(s{1}+1:end));
        phs = Descrip(Positions(2));
        s = strfind(phs,' ');
        obs(i,cc*2) = str2double(phs{1}(s{1}+1:end));
    end
end
% Make names
for ii = 1:length(const)
    const_names{2*ii-1} = [const{ii} '_amp'];
    const_names{2*ii} = [const{ii} '_phs'];
end
T1 = array2table(obs,'VariableNames',const_names);
T.Description = [];
% Combine
T = [T T1];

%% Only keep stations within our domain
opedat = m.op;
boudat = m.bd;
VX     = m.p;
%idx = convhull(VX(:,1),VX(:,2));
%idx = boundary(x,y);

in = T.Lon > -88.5 & T.Lon < -55 & T.Lat > 9 & T.Lat < 24;

%in = inpoly([T.Lon,T.Lat],[VX(idx,1) VX(idx,2)]);

% poly = []; 
% boudat.nbvv = full(boudat.nbvv);
% for op = 1:opedat.nope
%     % Get first openboundary
%     nodes = opedat.nbdv(1:opedat.nvdll(op),op);
%     if exist('node_o','var')
%         if nodes(1) ~= node_o(end)
%             nodes = flipud(nodes);
%         end
%     end
%     poly = [poly;  VX(nodes,:)];
%     
%     % Add the next land boundary
%     for i = 1:boudat.nbou
%         if boudat.nbvv(1,i) == nodes(end) || ...
%            boudat.nbvv(boudat.nvell(i),i) == nodes(end)
%             if boudat.nbvv(1,i) == nodes(end)
%                 nodes = boudat.nbvv(1:boudat.nvell(i),i);
%             else
%                 nodes = flipud(boudat.nbvv(1:boudat.nvell(i),i));
%             end
%             node_o = nodes;
%             mainland(op) = i;
%             poly = [poly;  VX(nodes,:)];
%             break;
%         end
%     end
% end
% if ~ispolycw(poly(:,1),poly(:,2))
%   poly = flipud(poly);
% end
% poly = [poly; NaN NaN];
% % Get the remaining boundaries
% for i = 1:boudat.nbou
%     if any(i == mainland); continue; end
%     nodes = boudat.nbvv(1:boudat.nvell(i),i);
%     if ispolycw(VX(nodes,1),VX(nodes,2))
%         nodes = flipud(nodes);
%     end
%     poly = [poly;  VX(nodes,:); NaN NaN]; 
% end
% edges = Get_poly_edges(poly);
% in = inpoly([T.Lon,T.Lat],poly,edges);

T = T(in,:);

%% Delete duplicates
radius = 0.08/60; % only take one station within 1 min box
[~,IA] = uniquetol([T.Lon T.Lat],radius,'ByRows',true,...
                    'DataScale',1,'OutputAllIndices',true);
% Just initialise T_new randomly
T_new = T(1:length(IA),:);
delete = []; nn = 0;
for ii = 1:length(IA)
    T_temp = T(IA{ii},:);
    if strcmp(T_temp.Name{1},'Haing Gyi Kyun')
    end
    % We want to pick the source in heirachial order
    for s = sources
        exp_s = regexp(T_temp.Source,s{1});
        exp_s = find(~cellfun(@isempty,exp_s)); 
        if ~isempty(exp_s)
            if length(exp_s) > 1
            end
            T_new(ii,:) = T_temp(exp_s,:);
            break
        end
    end
    if isempty(exp_s)
       nn = nn + 1;
       delete(nn) = ii; 
    end
end
T_new(delete,:) = [];
% sorting
[~,IA] = sort(T_new.Name);
T_new = T_new(IA,:);
[~,IA] = sort(T_new.Source);
T_new = T_new(IA,:);


%% Write to the .csv
% Output table of the names and positions
writetable(T_new(IA,:),output)
