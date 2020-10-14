% Plots the tidal amplitude and phases of desired constituents from a list of fort.53.nc files
clearvars; clc; close all; 

addpath(genpath('~/MATLAB/m_map'))
addpath(genpath('~/MATLAB/OceanMesh2D'))

% list of fort 53 files
filenames = dir('*53.nc');

% plot projection
projection = 'Miller';

% the constituent(s) to look at
conn = {'M2'}; 

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
       m_tricontour(ele,[x,y],g_m,0:30:360,'k');
       ax = gca;
       ax.Color = [166 132 97]/255;
       m_grid()
       caxis([0 1.6])
       colormap(lansey(16))
       cb = colorbar;
       cb.Label.String = 'amplitude [m]';
       print([fort53(1:end-6) '_' conn{c} '_amp+phase.png'],'-dpng','-r300')
   end
end
