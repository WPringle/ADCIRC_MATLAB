%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tidal analysis of data at all stations in the Tide_locations_dates.csv  %
% and output the constituent amplitudes and phases in a .cvs table file   %
% and a .psxy file to produce a figgen scatter.                           %
%                                                                         %
% Requires: Utide functions, obtained from;                               %
% https://www.mathworks.com/matlabcentral/fileexchange/                   %
% 46523--utide--unified-tidal-analysis-and-prediction-functions           %
%                                                                         %
% NETCDF data obtained from: http://uhslc.soest.hawaii.edu/data/?rq       %
%                                                                         %
% William Pringle Dec 2016                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializations
clearvars; close all; clc;
const_n = {'Name','Country','Latitude','Longitude','Start','End',...
           'M2_amp','M2_phs','S2_amp','S2_phs',...
           'N2_amp','N2_phs','K2_amp','K2_phs',...
           'K1_amp','K1_phs','O1_amp','O1_phs',...
           'P1_amp','P1_phs','Q1_amp','Q1_phs'};
%
%fid = fopen('Alltides_HC_allyears.txt','w');
ys = 1980; ye = 2016;
year = datenum(['01-Jan-' num2str(ys)]); 
year_n = datenum(['01-Jan-' num2str(ye)]); 
year_ref = datenum('01-Jan-1700');
file_s = 'OS_UH-RQH';
file_R = '_20160323_R.nc';
file_D = '_20160323_D.nc';
M = importdata('Tide_locations_dates.csv');
L = length(M.textdata);
Lats     = str2double(M.textdata(2:L,6));
Long     = str2double(M.textdata(2:L,7));
Name     = M.textdata(2:L,1);
Location = M.textdata(2:L,4);
Country  = M.textdata(2:L,5);
Startdate = M.textdata(2:L,8);
Enddate   = M.textdata(2:L,9);
Version  = M.textdata(2:L,3);
Data     = M.data(:,1);
NetCDF   = M.data(:,3);
const    = {'M2','S2','N2','K2','K1','O1','P1','Q1'};
format   = cell(length(const)+1,1); 
for i = 1:length(format)
    format{i} = '%s ';
end
%fprintf(fid,[format{:} '\n'],'Location',const{:});
format   = cell(length(const)*2,1); 
for i = 1:length(const)*2
    format{i} = '%7.3f, ';
end

%% Doing the analysis
for s = 1:length(Name)
   station{s,1} = Location{s};
   station{s,2} = Country{s};
   station{s,3} = Lats(s);
   station{s,4} = Long(s);
   station{s,5} = Startdate{s};
   station{s,6} = Enddate{s};
    if NetCDF(s) == 1
        FileName = [file_s Name{s} upper(Version{s}) file_R];
        if ~exist(FileName, 'file')
            FileName = [file_s Name{s} upper(Version{s}) file_D];
        end
        ncid = netcdf.open(FileName,'NC_NOWRITE');
        gattname = netcdf.inqAttName(ncid,netcdf.getConstant('NC_GLOBAL'),7);
        % No need to change timezone - keep in GMT
%         if strcmp(gattname,'time_zone')
%             % We need to add on the timezone to get into 
%             % local time to compare with Austides
%             tz = ncreadatt(FileName,'/','time_zone'); tzl = length(tz);
%             dt1 = str2double(tz(5)); dt2 = str2double(tz(7)); dt = dt1 + dt2/6;
        time = ncread(FileName,'time');
%            time = time + dt/24; % add on time zone in days
        zeta = ncread(FileName,'sea_surface_height_above_reference_level');
        zeta = squeeze(zeta);
        % Get only the desired years' data
        %I = find(time < (year - year_ref) | time >= (year_n - year_ref));
        %time(I) = []; zeta(I) = []; 
        time = time + year_ref;
    elseif Data(s) == 1
        FileName = ['h' Name{s} Version{s} '/i' Name{s} Version{s}];
        ne = 0;
        clear time zeta
        for y = ys:ye
            num = num2str(y); num = num(3:4);
            FileNameN = [FileName num '.dat'];
            if exist(FileNameN, 'file')
                D = dlmread(FileNameN,'',1,2);
                timeN = D(:,1); timeN = num2str(timeN);
                timey = zeros(length(timeN)*12,1); 
                zetay = zeros(size(timey));
                for j = 1:length(timeN)
                    timey(12*(j-1)+1) = datenum(str2double(timeN(j,1:4)),...
                            str2double(timeN(j,5:6)),str2double(timeN(j,7:8)));
                    timey(12*(j-1)+1) = timey(12*(j-1)+1) + (str2double(timeN(j,9))-1)/2;
                    zetay(12*(j-1)+1:12*(j-1)+12) = D(j,2:end);
                    for k = 2:12
                        timey(12*(j-1)+k) = timey(12*(j-1)+1) + (k-1)/24;
                    end
                end
                ns = ne + 1;
                ne = ne + length(timey);
                time(ns:ne) = timey;
                zeta(ns:ne) = zetay;
                zeta(zeta == 9999) = NaN;
            end
        end
    end
    entry = zeros(1,16);
    if length(find(~isnan(zeta))) >= 2190 % three months
        coef = ut_solv ( time, zeta, [], Lats(s), const,...
                        'OrderCnstit',const,'RunTimeDisp','nnn');
        % Reorder coefficient and change amplitude to metres
        for k = 1:length(const)
            for p = 1:length(coef.name)
                if coef.name{p} == const{k}
                    station{s,6+2*k-1} = coef.A(p)/1000;
                    station{s,6+2*k} = coef.g(p);
                    entry(2*k-1) = coef.A(p)/1000;
                    entry(2*k)   = coef.g(p);
                    break
                end
            end
            if s == 1
                fida{k} = fopen([const{k} '_amp.63.nc.psxy'],'w');
                fidp{k} = fopen([const{k} '_phs.63.nc.psxy'],'w');
            end
            fprintf(fida{k},'%9.3f %9.3f %8.3f %s \n',...
                    Long(s),Lats(s),entry(2*k-1),'5p');
            fprintf(fidp{k},'%9.3f %9.3f %8.3f %s \n',...
                    Long(s),Lats(s),entry(2*k),'5p');
         end
%         fprintf(fid,['%s ,' format{:} '\n'],Location{s},entry);
    end
end

%% Printing into table
T = cell2table(station,'VariableNames',const_n);
writetable(T,'IndianOcean_HC_allyears.csv')