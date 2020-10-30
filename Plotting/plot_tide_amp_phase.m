% Plots the tidal amplitude and phases of desired constituents from a list of fort.53.nc files
clearvars; clc; close all; 

%% User sets up their inputs and parameters here
% adding your paths for m_map, OM2D and data etc
addpath(genpath('~/MATLAB/m_map'))
addpath(genpath('~/MATLAB/OceanMesh2D'))

% list of fort 53 files
filenames = dir('*53.nc');
% output folder (make empty if save to current one)
outdir = 'figs/';

% % some parameters
% plot projection
projection = 'Miller';
% the constituent(s) to look at
conn = {'M2'}; 
% color range for amp
crange = [0 1.6]; % [m]
% number of discrete colors in colormap 
ncolor = 16; 
% contour intervals
cint = 0:30:360; % every 30 deg
% background color (this is a sort of brown)       
bgc = [166 132 97]/255;

%% Computation starts here
% make output directory if doesn't exist
if ~isempty(outdir) && ~isfolder(outdir)
   mkdir(outdir)
end
% Loop over all the fort.53.nc ADCIRC files
for ff = 1:length(filenames)
   fort53 = filenames(ff).name; 
   disp(fort53)
   x = ncread(fort53,'x');
   y = ncread(fort53,'y');
   nn = length(x);
   ele = ncread(fort53,'element')';
   xt = [x(ele(:,1)) x(ele(:,2)) x(ele(:,3),1) x(ele(:,1),1)];
   dxt = diff(xt,[],2);
   % for global mesh
   ele(abs(dxt(:,1)) > 180 | abs(dxt(:,2)) > 180 | ...
       abs(dxt(:,3)) > 180,:) = [];
   const = string(strtrim(ncread(fort53,'const')'));
   for c = 1:length(conn)
       disp(conn{c}) 
       ii = find(strcmp(const,conn{c}))          
  
       % Get modeled amp and phase of required constituent
       a_m = ncread(fort53,'amp',[ii 1],[1 nn]);
       g_m = ncread(fort53,'phs',[ii 1],[1 nn]);
                
       figure;
       m_proj(projection,'long',[min(x) max(x)],'lat',[min(y) max(y)])
       m_trisurf(ele,x,y,a_m);
       m_tricontour(ele,[x,y],g_m,cint,'k');
       ax = gca;
       ax.Color = bgc;
       m_grid()
       caxis(crange)
       colormap(lansey(ncolor))
       cb = colorbar;
       cb.Label.String = 'amplitude [m]';
       print([outdir filenames(ff).name(1:end-6) '_' conn{c} '_amp+phase.png'],'-dpng','-r300')
   end
end
