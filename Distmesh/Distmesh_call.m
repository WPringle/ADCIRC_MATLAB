% Distmesh calling script
clearvars; close all; 

%% Setting parth and compling InPolygon if required
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB')
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB/distmesh1.2')
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB/m_map')
%path(path,'/afs/crc.nd.edu/user/w/wpringle/MATLAB/InPolygon-MEX')
%mex '/afs/crc.nd.edu/user/w/wpringle/MATLAB/InPolygon-MEX/InPolygon.c' 

%% Set up parallel pool
Proc_num = 4;
if isempty(gcp)
    parpool('local',Proc_num);
end

%% Set parameters
mapfile = 'IDIOMS_v7.map';                           
bathyfile = ['E:\Global_Data\SRTM30_PLUS_w_Abyssal_Hills\' ...
            'bathy_SSG_1_120_GLOBAL_landmask.nc'];
edgelength = 50e3;    % min edgelength in meters
dist_param = 0.1;      % Distance paramater (percent that the edgelength 
                       % should change with distance, set zero to ignore)
wl_param = 15;         % parameter in wavelength function (set zero to ignore)
slope_param = 15;      % parameter in slope function (set zero to ignore)
itmax       = 10;      % Maximum number of iterations allowed in distmesh
plot_on = 1;           % Plot? (Yes = 1, No = 0)
outputname  = 'IDIOMS_out'; % Output .mat name;
 
%% Call distmesh
tic
[p,t] = General_distmesh(mapfile,bathyfile,edgelength,...
                         dist_param,wl_param,slope_param,itmax,plot_on);
toc

%% Close parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

%% Save outputs to .mat file
save([outputname '_' num2str(edgelength*1e-3) 'km.mat'],'p','t');