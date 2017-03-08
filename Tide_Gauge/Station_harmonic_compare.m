% Read fort 61 and get harmonic analysis of the time series
% Compare with the observed stuff
clearvars; clc; close all;

% Choose const.
const = {'M2','K1'};

% error type (percentage or error)
%type = '% error';
type = 'error';

alim = [0.25, 0.05];

% Read the data of the tide guages
tdirec = 'E:\Global_Data\Tidal Station Data\UHSLC Hourly Data\';

kmlS = kml2struct([tdirec 'Global_tide_gauge_database.kml']);

% Get the observed timeseries signal out
year_ref = ncreadatt('fort.61.nc','time','base_date');
year_ref = datenum(year_ref(1:19));
time = ncread('fort.61.nc','time');
time = time/(24*3600) + year_ref;
zeta = ncread('fort.61.nc','zeta');
name = ncread('fort.61.nc','station_name');
lon = ncread('fort.61.nc','x');
lat = ncread('fort.61.nc','y');

% Find the gappiness of each signal
I = []; %[9:length(lon)];
for i = 1:length(lon)
    L = length(find(~isnan(zeta(i,:))));
    % Ignore if data quite gappy
    if L/size(zeta,2) < 0.8
        I = [I; i];
    end
end
% Delete gappy time signal
name(:,I) = []; lon(I) = []; lat(I) = []; zeta(I,:) = [];

% Get only matching kmlS stuff
k = 0;
for c = const
    k = k + 1;
    figure(k);
    m_proj('Miller Cylindrical','lon',[21 132],'lat',[-75 31]);
    m_coast;
    hold on
end
c_o = zeros(length(lon),length(const));
c_n = zeros(length(lon),length(const));
for i = 1:length(kmlS)
    [I,dist] = knnsearch([lon,lat],[kmlS(i).Lon kmlS(i).Lat]);
    if dist < 1d-5
%         plot(kmlS(i).Lon,kmlS(i).Lat,'x')
%         hold on
%         plot(lon(I),lat(I),'o')

        % Do the harmonic analysis of numeric data
        coef = ut_solv ( time, zeta(I,:)', [], lat(I), const,...
                         'OrderCnstit',const,'RunTimeDisp','nnn');        

        % plot errors
        k = 0;
        for c = const
            k = k + 1;
            figure(k);
            
            % The observed data
            for j = 1:length(kmlS(i).Description)
                if strcmp(kmlS(i).Description{j}(1:2),c{1})
                    obs = str2num(kmlS(i).Description{j}(8:end));
                    break;
                end
            end   
            if isempty(obs)
                continue;
            end
             % For getting the phase for calculating rms
            phs_o = str2num(kmlS(i).Description{j+1}(10:end));
            c_o(I,k) = obs*exp(1i*pi*phs_o/180);
            for cc = 1:length(const)
                if strcmp(coef.name(cc),c{1})
                    num = coef.A(cc);
                    phs_n = coef.g(cc);
                    break;
                end
            end
            c_n(I,k) = num*exp(1i*pi*phs_n/180);
            if strcmp(type,'% error')
               error = (num - obs)/obs;
            elseif strcmp(type,'error')
               error = num - obs;
            end  
            [X,Y] = m_ll2xy(lon(I),lat(I));
            scatter(X,Y,50,error,'filled','MarkerEdgeColor','k')  
        end
    end
end
k = 0;
for c = const
    k = k + 1;
    figure(k);
    colormap(redblue)
    if strcmp(type,'% error')
        caxis([-1 1])
    elseif strcmp(type,'error')
        caxis([-alim(k),alim(k)])
    end
    colorbar
    title([c{1} ' amplitude ' type])
    rms_c = rms(c_n(:,k) - c_o(:,k));
    m_text(80,40,['RMS_{cmp} = ' num2str(rms_c*100) ' cm'])
end