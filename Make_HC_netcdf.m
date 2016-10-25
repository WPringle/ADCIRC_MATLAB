%% Open and plot output nc files
clear; clc;
fileID = fopen('fort.53');
T = fscanf(fileID,'%d',1);
C = textscan(fileID,'%f %f %f %s \n',T); C = C{4};
fclose(fileID);
fort53 = dlmread('fort.53','',T+1,0); F = length(fort53);
HC_amp = zeros(fort53(1,1),1); HC_phs = zeros(fort53(1,1),1);
for k = 1:T
    HC_amp = fort53(2+k:T+1:F,1);
    HC_phs = fort53(2+k:T+1:F,2);
    
    filename = [C{k} '_amp.63.nc'];
    copyfile('maxele.63.nc',[C{k} '_amp.63.nc'])
    ncid = netcdf.open(filename,'NC_WRITE');
    netcdf.putVar(ncid,16,HC_amp);
    netcdf.close(ncid)
    
    filename = [C{k} '_phs.63.nc'];
    copyfile('maxele.63.nc',[C{k} '_phs.63.nc'])
    ncid = netcdf.open(filename,'NC_WRITE');
    netcdf.putVar(ncid,16,HC_phs);
    netcdf.close(ncid)
end
% filename = 'fort_HC_amp.63.nc';
% ncid = netcdf.open(filename,'NC_WRITE');
% netcdf.putVar(ncid,16,[0 0],[fort53(1,1) T],HC_amp);
% netcdf.close(ncid)
% filename = 'fort_HC_phs.63.nc';
% ncid = netcdf.open(filename,'NC_WRITE');
% netcdf.putVar(ncid,16,[0 0],[fort53(1,1) T],HC_phs);
% netcdf.close(ncid)