function [coefS,coefL] = two_step_harm_analysis(time,zeta,lat,SNRmin)
% [coefS,coefL] = two_step_harm_analysis(time,zeta,SNRmin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tidal analysis of data of elevation time series using two-step process  %
%                                                                         %
% Requires: Utide functions, obtained from;                               %
% https://www.mathworks.com/matlabcentral/fileexchange/                   %
% 46523--utide--unified-tidal-analysis-and-prediction-functions           %
%                                                                         %
% Hourly time series, e.g. GESLA-2 website: http://gesla.org/             %
%                                                                         %
% William Pringle Nov 2017                                                %
%                                                                         %
% Inputs                                                                  %
%       time  : vector of time in datenum format                          %
%       zeta  : vector of elevations                                      %
%       lat   : latitude of the station                                   %
%       SNRmin: minimum signal-to-noise ratio for diagnostics             %
%                                                                         %
% Outputs                                                                 %
%       coefS : structure output from Utide for the monthly moving mean   %
%               averaged time series (long-term should not be valid)      %
%       coefL : structure output from Utide of the residual time series   %
%               after subtracting reconstructured time series from        %
%               coefS result to get long-term constituents                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% Here we try to detect sudden changes in constituent values between
%% years to trim time-series if needed
% Getting the matrix of times and elevations for each year
TE = datetime(time(end), 'ConvertFrom', 'datenum') ;
TS = datetime(time(1), 'ConvertFrom', 'datenum') ;
for i = 1:floor(years(TE-TS))
    tn = time(time <= datenum(TE-years(i-1)) & ...
              time > datenum(TE-years(i)));
    zn = zeta(time <= datenum(TE-years(i-1)) & ...
              time > datenum(TE-years(i)));         
    TM(1:length(tn),i) = tn;       
    ZM(1:length(tn),i) = zn;    
end
ZM(ZM == 0) = NaN; TM(TM == 0) = NaN; 

% Just getting the four main constituents for each year
coefM = ut_solv ( TM, ZM, [], M.OBS_LAT(ii),{'M2','S2','K1','O1'},...
                 'NoTrend');  
             
% Computing the semidiurnal/diurnal ratio
TNum = (coefM.A(1,:) + coefM.A(2,:))./(coefM.A(3,:) + coefM.A(4,:));
if mean(TNum) > 1
    TNum = 1./TNum;
end

% Detecting years with significant change
n = 0; I = [];
for iii = 1:length(TNum)-1
    TM1 = abs(mean(TNum(1:iii)) - mean(TNum(iii+1:end)));
    if abs(TM1) > 0.05
        n = n + 1;
        % long term 5% change
        I(n) = iii;
    end
end
% Shortening time series with significant change
if ~isempty(I)
    zeta = zeta(time > datenum(TE-years(I(1))));
    time = time(time > datenum(TE-years(I(1))));
end

%% Starting first-step
% Remove MSL with 30 day running mean so long-term 
% constituents become unimportant
zeta_orig = zeta;
zeta = zeta - movmean(zeta,30,'omitnan','SamplePoints',time); 

% Perform the analysis for all constituents 
% (will not give reliable estimates for the long-term constituents)
try 
    coefS = ut_solv ( time, zeta, [], lat, 'auto',...
                     'RunTimeDisp','nnn','NoTrend',...
                     'DiagnMinSNR', SNRmin);  
catch
    coefS = ut_solv ( time, zeta, [], lat, 'auto',...
                     'RunTimeDisp','nnn','NoTrend','LinCI',...
                     'DiagnMinSNR', SNRmin);  
end

%% Starting second-step - Long term
% Zero out long term and other low SNR constituents for reconstruction
coefR = coefS;
coefR.A(coef.diagn.SNR < SNRmin) = 0; 
coefR.A_ci(coef.diagn.SNR < SNRmin) = 0;
coefR.g(coef.diagn.SNR < SNRmin) = 0; 
coefR.g_ci(coef.diagn.SNR < SNRmin) = 0;
week = 1/(24*7); % One week frequency
coefR.A(coef.aux.frq < week) = 0; coefR.A_ci(coef.diagn.SNR < week) = 0;
coefR.g(coef.aux.frq < week) = 0; coefR.g_ci(coef.diagn.SNR < week) = 0;
% Reconstruct time series
sl_fit = ut_reconstr(time,coefR);
% Getting residual
residual = zeta_orig - sl_fit;

% Doing the analysis
try 
    coefL = ut_solv ( time, residual, [], lat,...
                  coef.name(coef.aux.frq < week),...
                  'RunTimeDisp','nnn','DiagnMinSNR', SNRmin);   
catch
    coefL = ut_solv ( time, residual, [], lat,...
                  coef.name(coef.aux.frq < week),...
                  'RunTimeDisp','nnn','LinCI','DiagnMinSNR', SNRmin);    
end
% EOF  
end