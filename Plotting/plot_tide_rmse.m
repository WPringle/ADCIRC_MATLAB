% Plots the tidal constituent rmse versus TPXO solutions
% from a list of fort.53.nc files
clearvars; clc; close all; 

%% User sets up their inputs and parameters here
% adding your paths for m_map, OM2D and data etc
addpath(genpath('~/MATLAB/m_map'))
addpath(genpath('~/MATLAB/OceanMesh2D'))

% list of fort 53 files
filenames = dir('*53.nc');
% output folder (make empty if save to current one)
outdir = 'figs/';

% TPXO9 reference file 
% (keep the ** wildcards and the script will replace with the constituent name)
tpxo9 = 'h_**_tpxo9_atlas_30.nc';

% % some parameters
% plot projection
projection = 'Miller';
% the constituent(s) to look at
conn = {'M2'}; 
% color range for rmse
crange = [0 0.25]; % [m]
% number of discrete colors in colormap 
ncolor = 20; 
% background color (this is a sort of brown)       
bgc = [166 132 97]/255;

%% Computation starts here
% make output directory if doesn't exist
if ~isempty(outdir) && ~isfolder(outdir)
   mkdir(outdir)
end
% Read the tpxo file attributes
istar = strfind(tpxo9,'*'); 
if ~isempty(istar); tpxo9(istar:istar+1) = 'm2'; end
lon = ncread(tpxo9,'lon_z');
lat = ncread(tpxo9,'lat_z');
if size(lon,2) == 1
    [lon ,lat] = meshgrid(lon,lat);
end
ctpx   = string(upper(strtrim(ncread(tpxo9,'con')')));
if ~isempty(istar); tpxo9(istar:istar+1) = '**'; end

% Loop over all the fort.53.nc ADCIRC files
for ff = 1:length(filenames)
   fort53 = [filenames(ff).folder '/' filenames(ff).name]; 
   disp(fort53)
   % Read the adcirc fort.53.nc attributes
   x = ncread(fort53,'x');
   y = ncread(fort53,'y');
   nn = length(x);
   ele = ncread(fort53,'element')';
   xt = [x(ele(:,1)) x(ele(:,2)) x(ele(:,3)) x(ele(:,1))];
   dxt = diff(xt,[],2);
   % get element mid-points: making sure we wrap around
   %xe = x(ele(:,1:3)); 
   %I = abs(dxt(:,1)) > 180 & abs(dxt(:,3)) > 180;
   %xe(I,1) = xe(I,1) + sign(dxt(I,1))*360;
   %I = abs(dxt(:,1)) > 180 & abs(dxt(:,2)) > 180;
   %xe(I,2) = xe(I,2) + sign(dxt(I,2))*360;
   %I = abs(dxt(:,2)) > 180 & abs(dxt(:,3)) > 180;
   %xe(I,3) = xe(I,3) + sign(dxt(I,3))*360;
   %xc = mean(xe,2);
   % put lon in 0-360 format for comparing to tpxo 
   %xc(xc < 0) = xc(xc < 0) + 360;
   xx = x; xx(xx < 0) = xx(xx < 0) + 360;
   %ye = y(ele(:,1:3)); yc = mean(ye,2);
   %clear xe ye xt dxt
   % for global mesh plotting
   ele(abs(dxt(:,1)) > 180 | abs(dxt(:,2)) > 180 | ...
       abs(dxt(:,3)) > 180,:) = [];
   const = string(strtrim(ncread(fort53,'const')'));

   for c = 1:length(conn)
       disp(conn{c}) 
       ii = find(strcmp(const,conn{c}));     
  
       % Get modeled amp and phase of required constituent
       a_m = ncread(fort53,'amp',[ii 1],[1 nn]);
       g_m = ncread(fort53,'phs',[ii 1],[1 nn]);
       a_m = a_m';
       g_m(g_m > 180) = g_m(g_m > 180) - 360;
       g_m = deg2rad(g_m');
       %% interpolate to the element mid-point
       %Z = a_m.*exp(1i*g_m);
       %F = scatteredInterpolant(x,y,Z);
       %Z = F(xc,yc);
       %a_m = abs(Z); g_m = angle(Z);
                
       % Get tpxo data of require constituent  
       istar = strfind(tpxo9,'*'); 
       if ~isempty(istar); 
          tpxo9(istar:istar+1) = lower(conn{c})
          Re = double(ncread(tpxo9,'hRe'))/1000;
          Im = double(ncread(tpxo9,'hIm'))/1000;  
          tpxo9(istar:istar+1) = '**';
       else
          jj = find(contains(ctpx,conn{c}));
          Re = ncread(tpxo9,'hRe',[1 1 jj],[size(lon) 1]);
          Im = ncread(tpxo9,'hIm',[1 1 jj],[size(lon) 1]);
       end
       % Interpolate to the element midpoint
       Z = Re - Im*1i;
       F = griddedInterpolant(lon',lat',transpose(Z));
       %BZ = F(xc,yc);
       BZ = F(xx,y);
       a_o = abs(BZ);  
       g_o = angle(BZ);
            
       % Compute the tidal rms difference
       RMS = sqrt(0.5*(a_o.^2 + a_m.^2 - 2*a_o.*a_m.*cos(g_o - g_m)));

       figure;
       m_proj(projection,'long',[min(x) max(x)],'lat',[min(y) max(y)])
       m_trisurf(ele,x,y,real(RMS));
       ax = gca;
       ax.Color = bgc;
       m_grid()
       caxis(crange)
       cmocean('amp',ncolor)
       cb = colorbar;
       cb.Label.String = 'RMSE [m]';
       print([outdir '/' filenames(ff).name(1:end-6) '_' conn{c} '_rmse.png'],'-dpng','-r300')
   end
end
