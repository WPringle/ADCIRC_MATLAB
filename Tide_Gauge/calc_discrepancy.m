% Calc RMS/total discrepancy
clearvars; close all

type_cell = {'Pelagic','Shelf','coast'};
%type_cell = {'coast'};

markertype = {'^','sq','o'};

fort53 = 'fort.53.nc';

c_list = {'M2','S2','N2','K2','K1','O1','P1','Q1'};

const_cell = {'M2','K1','all'};

x = ncread(fort53,'x');
y = ncread(fort53,'y');
boundary = ncread(fort53,'nbvv');
pb = [x(boundary),y(boundary)];
mdl1 = KDTreeSearcher(pb);
mdl2 = KDTreeSearcher([x,y]);
x = x'; y = y';

figure(1); hold on

c_list2 = ncread(fort53,'const');

tt = 0;
for type = type_cell 
    tt = tt + 1;
    T_o = readtable('ECGC_StationList_NoRepeats.csv');
    if strcmp(type{1},'coast')
        I1 = strfind(T_o.Source,'Truth');
        I2 = strfind(T_o.Name,'_TP');
        I3 = strfind(T_o.Name,'_Shelf');
        I = find(cellfun(@isempty,I1) & cellfun(@isempty,I2) & cellfun(@isempty,I3)); 
    else
        I1 = strfind(T_o.Source,type{1});
        I2 = strfind(T_o.Name,type{1});
        I = find(~cellfun(@isempty,I1) | ~cellfun(@isempty,I2));    
    end

    T_o = T_o(I,:);
    
    idx1 = knnsearch(mdl2,[T_o.Lon T_o.Lat],'k',12); 
    [~,distance] = knnsearch(mdl1,[T_o.Lon T_o.Lat]); 
    distance = distance*111e3;
    cc = 0;
    for const = const_cell
        cc = cc + 1;
        pos = 0; 
        RMS = zeros(length(I),1); V = zeros(length(I),1);
        for c = c_list
            pos = pos + 1;
            if ~strcmp(const{1},'all') && ~strcmp(c{1},const{1})
                continue
            end
            a_0 = table2array(T_o(:,4+2*pos-1));
            g_0 = table2array(T_o(:,4+2*pos));

            pos2 = 0;
            for ccc = c_list2
                pos2 = pos2 + 1;
                if strcmp(strtrim(ccc'),c{1})
                    break
                end
            end    

            % Using fort.53
            amp1 = ncread(fort53,'amp',[pos2 1],[1 length(x)]);
            phs1 = ncread(fort53,'phs',[pos2 1],[1 length(x)]);

            % get mode of phase (assume to be the null value - and make NaN)
            Pmode = mode(phs1(:));
            amp1(abs(phs1 - Pmode) < 1e-4) = NaN;
            phs1(abs(phs1 - Pmode) < 1e-4) = NaN;

            % These are nodes to make the interpolation
            [a_m, g_m] = sta_interp(x(idx1),y(idx1),amp1(idx1),...
                                      phs1(idx1),T_o.Lon,T_o.Lat);
                                  
            % Do not add to sum
            a_m(isnan(a_0)) = 0; g_m(isnan(a_0)) = 0;
            g_0(isnan(a_0)) = 0; a_0(isnan(a_0)) = 0; 
            g_0(isnan(a_m)) = 0; a_0(isnan(a_m)) = 0; 
            g_m(isnan(a_m)) = 0; a_m(isnan(a_m)) = 0; 
            
            %RMS = RMS + (a_0.^2 + a_m.^2 - 2*a_0.*a_m.*cosd(g_0 - g_m));
            RMS = RMS + (a_0.^2 + a_m.^2 - 2*a_0.*a_m.*cosd(g_0 - g_m));
            V   = V   + a_0.^2;
        end  
        K = find(RMS == 0);
        latt = T_o.Lat; lont = T_o.Lon;
        latt(K) = []; lont(K) = []; 
        RMS(K) = []; V (K) = [];
        %RSS = sqrt(0.5*sum(RMS)/length(RMS))*100
        RMS = sqrt(0.5*RMS); 
        V = sqrt(0.5*V); 
        
        D(tt,cc) = 100*mean(RMS);
        SD(tt,cc) = 100*std(RMS);
        %disp(D*100)

        RD(tt,cc) = 100*mean(RMS./V);
        SRD(tt,cc) = 100*std(RMS./V);
        
        D_cell{tt,cc} = RMS;
        V_cell{tt,cc} = V;
        Lat_cell{tt,cc} = latt;
        Lon_cell{tt,cc} = lont;
        
        if cc == 1
            figure(1);
            a_m(K) = []; a_0(K) = [];
            %s = scatter(Lon_cell{tt,1},Lat_cell{tt,1},[],D_cell{tt,1},'filled');
            s = scatter(Lon_cell{tt,cc},Lat_cell{tt,cc},[],a_m-a_0,'filled');
            s.Marker = markertype{tt};
        end
    end
end
plot_google_map('MapType','satellite')
%save('NI2/Discrepancy_results_NI2.mat','type_cell','const_cell',...
%                               'D_cell','V_cell','Lat_cell','Lon_cell',...
%                               'Dc_cell','Vc_cell','D','SD','RD','SRD',...
%                               'Dcoast','SDcoast','RDcoast','SRDcoast')
% T_o(isnan(RMS1) | RMS1 == 0,:) = [];
% V(isnan(RMS1) | RMS1 == 0) = [];
% %RMS(isnan(RMS1) | RMS1 == 0) = [];
% RMS1(isnan(RMS1) | RMS1 == 0) = [];
% 
% %RMS = sqrt(0.5*RMS); 
% RSS = sqrt(0.5*mean(RMS1));
% RMS1 = sqrt(0.5*RMS1); 
% V = sqrt(0.5*V); 
% 
% figure; hold on
% scatter(T_o.Lon,T_o.Lat,[],RMS1./V,'filled')
% caxis([0 0.5])
% 
% figure; hold on
% scatter(T_o.Lon,T_o.Lat,[],RMS1,'filled')
% caxis([0 0.25])
% 
% disp(RSS*100)
% 
% D = mean(RMS1);
% disp(D*100)
% 
% RD = mean(RMS1./V);
% disp(RD*100)



