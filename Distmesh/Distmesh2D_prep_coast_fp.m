function [p,t] = Distmesh2D_prep_coast_fp(meshfile,contourfile,bathyfile,...
                                 bbox,min_el,max_el,dist_param,wl_param,...
                                 slope_param,bounds,minL,itmax,...
                                 plot_on,ini_p,fix_p,nscreen)
% Function for calling grid generator for any general region in the world
% Inputs:  meshfile   : filename and path of mesh with mainland boundary
%          contourfile: filename(s) of shape file(s) (must be a cell)
%          bathyfile  : filename of the bathymetric dataset (the version here
%                       is for a global dataset in netcdf format)
%          bbox       : the bounding box to mesh within in following format:
%                       [minlon maxlon;
%                        minlat maxlat];
%          min_el     : The minimum (and initial) edgelength of the mesh
%          max_el     : The maximum edgelength of the mesh
%          dist_param : The fraction that the edgelength should change
%                       with distance. Set 0 to ignore
%          wl_param   : Parameter in wavelength function, (set 0 to ignore)
%                       edgelength = T/wl_param*sqrt(gH), where T is the
%                       period of the M2 wave and H is the bathymetric depth
%          slope_param: Parameter in slope function, (set 0 to ignore)
%                       edgelength = 2*pi/slope_param*H/abs(grad(H)),
%                       where H is the bathymetric depth
%          bounds     : Where meshfile is present (to mesh floodplain) 
%                       bounds(1) is the upper bound of the topo to mesh
%                       within and bounds(2) is the max dist from coastline
%                       to mesh within 
%          minL       : the minimum length of polygon points (island size) 
%                       to allow (usually 5-7)
%          itmax      : Maximum number of iterations in distmesh function
%          plot_on    : choose whether to plot map and mesh during calc or
%                       not (Yes = 1, No = 0)
%          ini_p      : optional initial distribution of p
%          fix_p      : optional fixed node points
%          n_screen   : Frequency of iterations to save a temporary mesh 
%                       and display the mesh quality
%
% Outputs: p          : the node positions in lon,lat
%          t          : the triangulation of the mesh
%
% Requires            : readfort14 (to read an ADCIRC mesh - floodplain only)
%                     : extdom_edges (to get edges of mesh - floodplain only)
%                     : extdom_polygon
%
% Initial development by William Pringle, Damrongsak Wirasaet, 2016-Feb 2017
% Combining floodplain+coastal meshing, parallel improvements, and 
% polygonal selection tool for floodplain meshing by Keith Roberts June 2017
%
% Reference for wavelength and slope function:
% Lyard, F., Lefevre, F., Letellier, T., & Francis, O. (2006).
% Modelling the global ocean tides: modern insights from FES2004.
% Ocean Dynamics, 56(5�6), 394�415. http://doi.org/10.1007/s10236-006-0086-x
%
%% Read the mesh if it is present (this assumes we are meshing the floodplain.)
if ~isempty(meshfile)
    tic
    disp('Reading the mesh...');
    [ev,pv,~,opendat,boudat,~] = readfort14(meshfile);
    disp('Read the mesh.');
    toc
    
    tic
    disp('Identifying the external domain segments...');
    segment = []; polygon = [];
    etbv    = extdom_edges( ev, pv ) ; % the edges in the external domain
    iedbeg  = knnsearch(pv(etbv(1,:),:),pv(opendat.nbdv(1),:)); % find which edge to start
    polygon = extdom_polygon( etbv, pv, iedbeg, 1, 2 ) ; % find the outer polygon defining the domain
    disp('Finished identifying the external domain segements.');
    toc
    
    tic
    disp('Creating the land boundary segments...');
    % Make boudat into segment and polygons based on min criteria
    for i = 1:boudat.nbou
        nodes = boudat.nbvv(1:boudat.nvell(i),i);
        if boudat.nvell(i) > minL
            % Add to segment
            segment = [segment; NaN NaN; pv(nodes,:)];
            %else
            %    % Add to polygon
            if ~any(polygon(:,1) == pv(nodes(1),1))
                if ~ispolycw(pv(nodes,1),pv(nodes,2))
                    polygon = [polygon; NaN NaN; pv(nodes,:)];
                else
                    polygon = [polygon; NaN NaN; flipud(pv(nodes,:))];
                end
            end
        end
    end
    % make fixp = to the segment
    fix_p  = segment;
    fix_p  = fix_p(~any(isnan(fix_p),2),:);
    toc
else
    polygon_struct = [];
    polygon        = [];
    % Read polygon from shape file and
    tic
    disp('Reading in shapefile...')
    polygon_struct = Read_shapefile( contourfile, bbox, minL, ...
                                     min_el, 0, polygon );
    
    % Smooth the polygon using a 5-point moving average
    polygon_struct.outer = smooth_coastline(polygon_struct.outer,5,plot_on);
    polygon_struct.mainland = smooth_coastline(polygon_struct.mainland,5,plot_on);
    polygon_struct.inner = smooth_coastline(polygon_struct.inner,5,plot_on);
    toc 
    
    % Plot the polygons
    if plot_on >= 1
        disp('Check bounding box...OK?'); pause;
    end
end

%% Read DEM and build bathy interpolant for fd (for floodplain)
tic
disp('Reading the DEM...');
lon = double(ncread(bathyfile,'lon'));
lat = double(ncread(bathyfile,'lat'));
I = find(lon > bbox(1,1) & lon < bbox(1,2));
J = find(lat > bbox(2,1) & lat < bbox(2,2));
lon = lon(I); lat = lat(J);
% read only part of the DEM necessary
ncid = netcdf.open(bathyfile,'NC_NOWRITE');
bathy = netcdf.getVar(ncid,2,[I(1) J(1)],[length(I) length(J)]);
netcdf.close(ncid)
%
[lon_g,lat_g] = ndgrid(lon,lat);
% Save the depth interpolant
Fb = griddedInterpolant(lon_g,lat_g,bathy);
toc

%% Make the distance kdtree searchers and interpolants 
if ~isempty(meshfile) %% Floodplain mesh
    % For distance function used to calculate edgelength function
    disp('Building KD-Tree (mdl1) with land boundary segments...');
    mdl1 = KDTreeSearcher(segment);
    
    % For floodplain prompt to draw a smaller bounding polygon to mesh within
    figure;
    plot(lon_g(1:100:end,1:100:end),lat_g(1:100:end,1:100:end),'m.') % bathymetry data availability within bbox
    hold on; plot(segment(:,1),segment(:,2),'r'); % mainland boundary
    disp('Please draw the outer bounding polygon that partially encompasses the mainland boundary...');
    h = impoly();
    floodplain_polygon = h.getPosition;
    mdl0 = KDTreeSearcher([segment;floodplain_polygon]); % this is used for inpoly
    
else %% Coastal Mesh
    % For the distance function used to calculate edgelength function
    disp('Building KD-Tree (mdl1) with mainland and island segments...');
    mdl1 = KDTreeSearcher([polygon_struct.inner; polygon_struct.mainland]);

    % For distance function used to determine whether points and elements
    % are inside the region you want to mesh
    disp('Building KD-Tree (mdl0) with ocean boundary, mainland, and island segments...');
    mdl0 = KDTreeSearcher([polygon_struct.outer; polygon_struct.inner]); 
    
%     % Build distance from outer boundary interpolant 
%     % (used to speed up distance evaluations). Note the difference from
%     % edgefunction using instance 0 here
%     x_v = reshape(lon_g,[],1);
%     y_v = reshape(lat_g,[],1);
%     disp('Calculating the distance from the ocean boundary...')
%     dOB = fd( [x_v,y_v], 0 ) ; %dpoly([x_v,y_v],[polygon_struct.outer; NaN,NaN; polygon_struct.inner],mdl0);
%     dOB = reshape(dOB,length(I),[]);
%     FOB = griddedInterpolant(lon_g,lat_g,dOB,'linear');
%     
%     % plot the distance function
%     if plot_on == 2
%         figure; contourf(lon_g,lat_g,dOB);
%         clearvars x_v y_v dOB
%     end
end

%% Make edge function interpolant
disp('Building edge interpolant...');
tic
nn = 0; % Counter for edge function options
% Distance..
if dist_param > 0
    disp('   Building distance edge function...')
    nn = nn + 1;
    x_v = reshape(lon_g,[],1);
    y_v = reshape(lat_g,[],1);
    d = fd( [x_v,y_v], 1 ) ;
    % reshape back
    d = reshape(d,length(I),[]);
    % add into edge function
    hh(:,:,nn) = min_el - dist_param*d ;
    %figure; contourf(lon_g,lat_g,hh(:,:,1));
    clear x_v y_v d
    % Plot the edgelength function for distance
    if plot_on == 2
        figure;
        m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)])
        m_contourf(lon_g,lat_g,hh(:,:,1)); shading interp
        hold on; m_plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2),'k-','linewi',2); 
        cb = colorbar; ylabel(cb,'Degrees from land'); 
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);    
        title('Distance function');
        print('Distancefunction','-dpng','-r300')
    end
end
% Wavelength...
if wl_param > 0
    disp('   Building the wavelength edge function...');
    nn = nn + 1;
    g = 9.807;
    period = 12.42*3600; % M2 period in seconds
    % wavelength func
    hh(:,:,nn) = period*sqrt(g*abs(bathy))/wl_param;
    % Plotting the wavelength edgelength function
    if plot_on == 2
        figure;
        m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)]);
        m_pcolor(lon_g,lat_g,real(hh(:,:,2))); shading interp
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        title('Wavelength function'); colorbar; 
        %ylabel(cb,'Elements per M2 wavelength'); 
        print('Wavelengthfunction','-dpng','-r300')
        hh(:,:,nn) = real(hh(:,:,nn));
    end
end
% Slope...
if slope_param > 0
    disp('   Building the slope edge function...');
    nn = nn + 1;
    
    % Converting lon and lat to metres so we can evaluate the slope
    % Note: Parallelization of m_idist which is very slow
    lon2 = [lon_g(2:end,:); lon_g(end,:)+(lon(end)-lon(end-1))];
    lat2 = [lat_g(:,2:end) lat_g(:,end)+(lat(end)-lat(end-1))];
    tic
    % lon
    dx  = zeros(size(lon2));
    np1 = size(lon2,1) ;
    Pool = gcp('nocreate'); num_p = Pool.NumWorkers;
    for idx = 1:num_p
        ns1 = int64((idx-1)*np1/num_p)+1;
        ne1 = int64(idx*np1/num_p);
        f(idx) = parfeval(Pool,@m_idist,2,...
                          lon_g(ns1:ne1,1:end),lat_g(ns1:ne1,1:end),...
                          lon2(ns1:ne1,1:end),lat_g(ns1:ne1,1:end));
    end
    for idx = 1:num_p
        [idx_t, dx_t, ~] = fetchNext(f); % Get results into a cell array
        ns1 = int64((idx_t-1)*np1/num_p)+1;
        ne1 = int64(idx_t*np1/num_p);
        dx(ns1:ne1,1:end) = dx_t;
    end
    % lat
    dy  = zeros(size(lon2));
    np1 = size(lon2,1) ;
    for idx = 1:num_p
        ns1 = int64((idx-1)*np1/num_p)+1;
        ne1 = int64(idx*np1/num_p);
        f(idx) = parfeval(Pool,@m_idist,2,...
                          lon_g(ns1:ne1,1:end),lat_g(ns1:ne1,1:end),...
                          lon_g(ns1:ne1,1:end),lat2(ns1:ne1,1:end));
    end
    for idx = 1:num_p
        [idx_t, dy_t, ~] = fetchNext(f); % Get results into a cell array
        ns1 = int64((idx_t-1)*np1/num_p)+1;
        ne1 = int64(idx_t*np1/num_p);
        dy(ns1:ne1,1:end) = dy_t;
    end
    toc
    % Get the coordinate matrices for dx and dy
    dx = cumsum(dx,1); dx = cumsum(dx,2);
    dy = cumsum(dy,1); dy = cumsum(dy,2);
    
    % Calculate the gradient of bathy
    [b_y,b_x] = gradient(bathy,dy,dx);
    b_slope = sqrt(b_x.^2 + b_y.^2);
    
    % slope func
    hh(:,:,nn) = 2*pi*bathy./b_slope/slope_param;
    clear lon2 lat2 dx dy b_y b_x b_slope
end

% Get min of slope and wavelength, setting equal to hh_m
disp('Finalizing the edge function..');
if wl_param > 0 || slope_param > 0
    if dist_param > 0
        hh_m = min(hh(:,:,2:end),[],3);
    else
        hh_m = min(hh(:,:,1:end),[],3);
    end
    % We need to convert hh_m into degrees from meters
    % (we estimate from 45 deg azimuth)
    % Note: Parallelized as very slow
    np1 = size(hh_m,1) ;
    lat2 = zeros(size(hh_m));
    lon2 = zeros(size(hh_m));
    for idx = 1:num_p
        ns1 = int64((idx-1)*np1/num_p)+1;
        ne1 = int64(idx*np1/num_p);
        f(idx) = parfeval(Pool,@m_fdist,3,...
                          lon_g(ns1:ne1,1:end),lat_g(ns1:ne1,1:end),...
                          45,hh_m(ns1:ne1,1:end));
    end
    for idx = 1:num_p
        [idx_t, lon2_t, lat2_t,~] = fetchNext(f); % Get results into a cell array
        ns1 = int64((idx_t-1)*np1/num_p)+1;
        ne1 = int64(idx_t*np1/num_p);
        lon2(ns1:ne1,1:end) = lon2_t;
        lat2(ns1:ne1,1:end) = lat2_t;
    end
    % switch back to -180 to 180
    lon2(lon2 > 180) = lon2(lon2 > 180) - 360;
    % get the actual distance
    hh_m = sqrt((lon2 - lon_g).^2 + (lat2 - lat_g).^2);
    
    % Get min of all the criteria
    if dist_param > 0
        hh_m = min(hh_m,hh(:,:,1));
    end
else
    hh_m = squeeze(hh);
end
clear lon2 lat2 hh

% Making sure edgelength not < minimum or > maximum allowable
hh_m(hh_m < min_el) = min_el; 
hh_m(hh_m > max_el) = max_el; 
disp('Finalized edge function');

if plot_on == 2
    figure; contourf(lon_g,lat_g,hh_m); 
    title('Before smoothing'); caxis([0,10*min_el]);colorbar;
end
% Smoothing edgelength function using the Lipschitz smoother from mesh2d
disp('Triangulating the bathymetric grid...')
x_v = reshape(lon_g,[],1); y_v = reshape(lat_g,[],1); 
pdmy = [x_v,y_v]; clear x_v y_v; 
tdmy = delaunay(pdmy(:,1),pdmy(:,2)); 
% We limit slope to be dist_param on the diagonal
hfun = reshape(hh_m,[],1); dhdx = dist_param/sqrt(2); 
% Enforce lipschitz continuity 
% to vary the gradient no more than the dist_param
[hfun,flag] = limhfn2(pdmy,tdmy,hfun,dhdx); 
if ~flag
   disp('Warning: was not able to satisfy smoothing criteria = dist_param')
end
hh_m = reshape(hfun,length(I),[]);
clear hfun pv pdmy tdmy % free memory 
if plot_on == 2
    figure; contourf(lon_g,lat_g,hh_m); 
    title('After Lipshitz smooth'); caxis([0,10*min_el]);colorbar;
end

% Showing the min edgelength in meters
avgLat = mean(mean(lat_g));
minLength = abs(111*1000*cos(avgLat)*min_el);
disp(['Minimum edge length is ',num2str(minLength),' meters.']);
    
% Plotting the overall edgelength function
if plot_on >= 1
    figure;
    m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)]);
    m_pcolor(lon_g,lat_g,hh_m); shading interp;
    m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
    title('Overall edgelength function'); colorbar;
    print('Overall_edgelength_function.png','-dpng','-r300')
end

% Make the overall interpolant
disp('Building the finalized edge function interpolant...')
F = griddedInterpolant(lon_g,lat_g,hh_m,'linear');
disp('Built the finalized edge function interpolant');
toc
clear hh_m lon_g lat_g

%% Call the distmesh iteration routine 
disp('Entering DistMesh2d');
close all;
[p,t] = distmesh2d(@fd,@fh,min_el,bbox',ini_p,fix_p,itmax,...
                   plot_on,nscreen);
disp('Exiting Distmesh2d.m');

%% Remove portion of the mesh that exceeds your distance AND depth criteria
if ~isempty(meshfile)
    disp('Trimmming mesh...')
    [p,t] = trim_mesh(p,t);
end

%% Fixing up the mesh automatically
disp('Fixing mesh...');
[p,t] = Fix_bad_edges_and_mesh(p,t);
tq = gettrimeshquan( p, t);
mq_m = mean(tq.qm);
mq_l = prctile(tq.qm,0.1);
disp(['number of nodes is ' num2str(length(p))])
disp(['mean quality is ' num2str(mq_m)])
disp(['0.1 prctile quality ' num2str(mq_l)])   
if plot_on >= 1 
    simpplot(p,t);
end
%
return;
%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% SUB_FUNCTIONS FOR DISTANCE, EDGELENGTH and TRIMMING    %  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
    function [p,t] = trim_mesh(p,t)
        % Trims the mesh (for floodplain only)
        pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3; % Compute centroids
        dis = fd(pmid,2);
        t = t(dis < 0,:);   % Keep interior triangle
        [p,t] = fixmesh(p,t); % clean up disjoint vertices.
        if plot_on >= 1
            simpplot(p,t);
        end
    end

    function d = fd( p, instance )
        % The distance function
        % needs segment to get distance
        % smallpolygons to make sure nodes inside them are discarded
        % Fb and bounds to test whether depth in desired range
        % NOTE: Pass mdl here instead of segment
        if instance == 0 
            % used for kicking points out of the box.
            % mdl0 contains the ocean boundary
            if ~isempty(meshfile)
                d = dpoly_fp2(p,mdl0,polygon,floodplain_polygon,Fb,[]);
            else
                d = dpoly(p,[polygon_struct.outer;
                             polygon_struct.inner],mdl0);
            end
        elseif instance == 1 
            % used for the creation of the distance function (sans ocean bou).
            % mdl1 does not contain the ocean boundary
            if ~isempty(meshfile)
                d = dpoly_fp2(p,mdl1,polygon,floodplain_polygon,[],[]);
            else
                d = dpoly(p,[polygon_struct.outer;
                             polygon_struct.inner],mdl1);
            end
        else
            % only used for floodplain meshing
            d = dpoly_fp2(p,mdl1,polygon,floodplain_polygon,Fb,bounds);
        end
        return ;
    end

    function h = fh( p )
        % The edge length function, desides on the edgelength at this
        % location, based on distance, wavelength and slope functions
        
        % Just call the interpolant created at beginning
        h = F(p);
        
        return ;
    end

end
