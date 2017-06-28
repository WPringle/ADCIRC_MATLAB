function [p,t] = Distmesh2D_prep(meshfile,contourfile,bathyfile,bbox,...
                                 min_el,max_el,dist_param,wl_param,...
                                 slope_param,bounds,minL,itmax,...
                                 plot_on,ini_p,fix_p,nscreen,num_p)
% Function for calling grid generator for any general region in the world
% Inputs:  meshfile   : filename and path of mesh with mainland boundary
%          contourfile: filename(s) of shape file(s) (must be a cell)
%          bathyfile  : filename of the bathymetric dataset (the version here
%                       is for a global dataset in netcdf format)
%          bbox       : the bounding box to mesh within in following format:
%                       [minlon maxlon;
%                        minlat maxlat];
%          min_el     : The minimum (and initial) edgelength of the mesh
%          max_el     : The maximum edgelength of the mesh. 
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
%          nscreen     : Frequency of iterations to save a temporary mesh
%                       and display the mesh quality
%          num_p      : number of parallel processors required 
%                      (<=1 for serial)
%
% Outputs: p          : the node positions in lon,lat
%          t          : the triangulation of the mesh
%
% Requires            : readfort14 (to read an ADCIRC mesh - floodplain only)
%                     : extdom_edges (to get edges of mesh - floodplain only)
%                     : extdom_polygon
%
% V0 : Initial development by William Pringle, Damrongsak Wirasaet, 2016-Feb 2017
% V1:  Combining floodplain+coastal meshing,overall improvements+bugfixes, and
%      polygonal selection tool for floodplain meshing by Keith Roberts 2017-April-2017-June
%
% Reference for wavelength and slope function:
% Lyard, F., Lefevre, F., Letellier, T., & Francis, O. (2006).
% Modelling the global ocean tides: modern insights from FES2004.
% Ocean Dynamics, 56(5�6), 394�415. http://doi.org/10.1007/s10236-006-0086-x
%
%% Read the mesh if it is present (this assumes we are meshing the floodplain.)
if ~isempty(meshfile)
    % Flooplain ----------------------------------------------------------
    disp('Reading the mesh...');
    [ev,pv,~,opendat,boudat,~] = readfort14(meshfile);
    
    disp('Identifying the external domain segments...');
    segment = []; polygon = [];
    % We first extract the outermost polygon
    etbv    = extdom_edges( ev, pv ) ; % the edges in the external domain
    iedbeg  = knnsearch(pv(etbv(1,:),:),pv(opendat.nbdv(1),:)); % find which edge to start
    polygon = extdom_polygon( etbv, pv, iedbeg, 1, 2 ) ; % find the outer polygon defining the domain
    
    disp('Creating the land boundary segments...');
    % Make boudat into segment and polygons based on min criteria
    for i = 1:boudat.nbou
        nodes = boudat.nbvv(1:boudat.nvell(i),i);
        if boudat.nvell(i) > minL
            % Add to segment
            segment = [segment; NaN NaN; pv(nodes,:)];
            % Add to polygon if doesn't already exist in outer one
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
    fixp  = segment;
    fixp  = fixp(~any(isnan(fixp),2),:);
    toc
elseif ~isempty(contourfile)
    % Coastal ----------------------------------------------------------
    boubox = [bbox(1,1) bbox(2,1); bbox(1,1) bbox(2,2); ...
              bbox(1,2) bbox(2,2); bbox(1,2) bbox(2,1); ...
              bbox(1,1) bbox(2,1)];
    polygon_struct = [];
    polygon        = [];
    % Read polygon from shape file and
    disp('Reading in shapefile...')
    polygon_struct = Read_shapefile( contourfile, bbox, minL, ...
                                     min_el, 0, polygon );
    % Smooth the polygons and segments that exist using a 5-point moving average
    poly_count = 0; oc = 0; inr = 0; ml = 0; ocean_only = 0;
    if ~isempty(polygon_struct.outer) 
        oc = oc + 1; poly_count = poly_count + 1; 
        [polygon_struct.outer,closed] = smooth_coastline(polygon_struct.outer,5,0);
    end
    if ~isempty(polygon_struct.mainland)
        inr = inr + 1; poly_count = poly_count + 1;
        polygon_struct.mainland = smooth_coastline(polygon_struct.mainland,5,0);
    end
    if ~isempty(polygon_struct.inner) 
        ml = ml + 1; poly_count = poly_count + 1;
        polygon_struct.inner = smooth_coastline(polygon_struct.inner,5,0);
    end
    % handle the case when no segments are present are in the bbox.
    if oc == 1 && poly_count == 1
        ocean_only = 1;
        disp('INFO: Meshing the open ocean, no land or island segments present within bbox...');
    end
    % template to remove vertical/horizontal artifacts from shapefile
    %     polygon_struct.outer(abs(polygon_struct.outer(:,1)- 78.00000) < 0.0001,1) = -9999;
    %     polygon_struct.outer(abs(polygon_struct.outer(:,1)- 78.00000) < 0.0001,2) = -9999;
    %     polygon_struct.outer(polygon_struct.outer(:,1)==-9999,:) = [];
%     if plot_on >= 1 
%         m_proj('Mercator','long',[bbox(1,1)-25, bbox(1,2)+25],...
%                           'lat',[bbox(2,1)-25, bbox(2,2)+25])
%         m_plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2),...
%                'm-','linewi',2); hold on;
%         m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
%         m_gshhs_i('color','b');       % Coastline...
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
%     end
    % Cornerify: create artifical outer bounding box to turn segment into
    % polygon. Lets test for this automatically
    if ~closed && ~ocean_only
        [polygon_struct.outer,ierr] = cornerify(polygon_struct.outer,boubox,plot_on);
        if ierr; return; end
    end
    if plot_on >=1
        m_proj('Mercator','long',[bbox(1,1)-25, bbox(1,2)+25],...
                          'lat', [bbox(2,1)-25, bbox(2,2)+25])
        hold on;
        % Plot the polygons that exist
        if inr == 1
           m_plot(polygon_struct.inner(:,1),polygon_struct.inner(:,2),'m');
        end
        if ml == 1
           m_plot(polygon_struct.mainland(:,1),polygon_struct.mainland(:,2),'k-.')
        end
        m_plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2),'k-','linewi',2); 
        m_hatch(boubox(:,1),boubox(:,2),'single',30,5,'color','k'); % ...with hatching added.
        m_gshhs_i('color','b'); % Coastline for reference...
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        title('Press 1 to accept or 0 to abort...');
        ierr = abort; if ierr==1; return; end
        save pg polygon_struct
    end
else % you have already run in interactive mode just load the polygon_struct
    disp(['INFO: Neither the contourfile or the meshfile were provided, '... 
         'loading in the polygon_struct from a previous interactive session...']);
    load pg
end

%% Read DEM and build bathy interpolant for fd (for floodplain)
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
[lon_g,lat_g] = ndgrid(lon,lat);
% Save the depth interpolant
Fb = griddedInterpolant(lon_g,lat_g,bathy);

%% Make the distance kdtree searchers and interpolants
if ~isempty(meshfile) %% Floodplain mesh ---------------------------------
    % For distance function used to calculate edgelength function
    disp('Building KD-Tree (mdl1) with land boundary segments...');
    mdl1 = KDTreeSearcher(segment);
    
    % For floodplain prompt to draw a smaller bounding polygon to mesh within
    figure;
    d_x = ceil(length(I)/100); % ensure only size 100 x 100
    d_y = ceil(length(J)/100); % ensure only size 100 x 100
    plot(lon_g(1:d_x:end,1:d_y:end),...
         lat_g(1:d_x:end,1:d_y:end),'m.') % bathymetry data availability within bbox
    hold on; plot(segment(:,1),segment(:,2),'r'); % mainland boundary
    disp(['Please draw the outer bounding polygon that ' ...
          'partially encompasses the mainland boundary...']);
    h = impoly();
    floodplain_polygon = h.getPosition;
    mdl0 = KDTreeSearcher([segment; floodplain_polygon]); % this is used for inpoly
    
else %% Coastal Mesh ------------------------------------------------------
    % For the distance function used to calculate edgelength function
    if ~ocean_only
        disp('Building KD-Tree (mdl1) with mainland and island segments...');
        mdl1 = KDTreeSearcher([polygon_struct.inner; polygon_struct.mainland]);
    else
        mdl1 = [];
    end
    
    % For distance function used to determine whether points and elements
    % are inside the region you want to mesh
    disp('Building KD-Tree (mdl0) with ocean boundary, mainland, and island segments...');
    if ~ocean_only
        mdl0 = KDTreeSearcher([polygon_struct.outer; NaN NaN; ...
                               polygon_struct.inner]);
    else
        mdl0 = KDTreeSearcher(polygon_struct.outer);
    end
    dOB = [];
    if size(bathy,1)*size(bathy,2) > 1e4 % if larger than e.g. 100 x 100
        dxx = abs(lon(2) - lon(1));
        % Build distance from outer boundary interpolant
        % (used to speed up distance evaluations for big areas or hi res data). 
        x_v = reshape(lon_g,[],1); y_v = reshape(lat_g,[],1);
        disp('Calculating the distance from the ocean boundary...');
        dOB = fd( [x_v,y_v], 0 ) ;
        dOB = reshape(dOB,length(I),[]);
        FOB = griddedInterpolant(lon_g,lat_g,dOB,'linear');
    end
        
    % plot the distance function for debug puproses to identify problems
    % with polygon closure.
    if plot_on == 3
        figure; contourf(lon_g,lat_g,dOB);
        title('Distance from outer polygon.');
        clearvars x_v y_v dOB
    end
end

%% Make edge function interpolant
disp('Building edge interpolant...');
nn = 0; % Counter for edge function options
% Distance..
if dist_param > 0 && ~ocean_only
    disp('   Building distance edge function...')
    nn = nn + 1;
    x_v = reshape(lon_g,[],1);
    y_v = reshape(lat_g,[],1);
    d = fd( [x_v,y_v], 1 ) ;
    % reshape back
    d = reshape(d,length(I),[]);
    % add into edge function
    hh(:,:,nn) = min_el - dist_param*d ;
    clear x_v y_v d
    % Plot the edgelength function for distance
    if plot_on >= 2
        figure;
        m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)])
        m_contourf(lon_g,lat_g,hh(:,:,1)); shading interp
        hold on; m_plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2),'k-','linewi',2);
        cb = colorbar; ylabel(cb,'Degrees from land');
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        title('Distance function');
        ierr = abort; if ierr == 1; return; end
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
    if plot_on >= 2
        figure;
        m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)]);
        m_pcolor(lon_g,lat_g,real(hh(:,:,nn))/1000); shading interp
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        title('Wavelength of the M2 edge function'); cb=colorbar;ylabel(cb,'km'); 
        ierr = abort; if ierr == 1; return; end
        print('Wavelengthfunction','-dpng','-r300')
    end
    hh(:,:,nn) = real(hh(:,:,nn));
end
% Slope...
if slope_param > 0
    disp('   Building the slope edge function...');
    nn = nn + 1;
    
    % Converting lon and lat to metres so we can evaluate the slope
    lon2 = [lon_g(2:end,:); lon_g(end,:)+(lon(end)-lon(end-1))];
    lat2 = [lat_g(:,2:end) lat_g(:,end)+(lat(end)-lat(end-1))];
    % lon
    if num_p > 1
        dx = m_idist_par(lon_g,lat_g,lon2,lat_g,num_p);
    else
        dx = m_idist(lon_g,lat_g,lon2,lat_g);
    end
    % lat
    if num_p > 1
        dy = m_idist_par(lon_g,lat_g,lon_g,lat2,num_p);
    else
        dy = m_idist(lon_g,lat_g,lon_g,lat2);
    end
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
if wl_param > 0 || slope_param > 0
    if dist_param > 0 && ~ocean_only
        hh_m = min(hh(:,:,2:end),[],3);
    else
        hh_m = min(hh(:,:,1:end),[],3);
    end
    % We need to convert hh_m into degrees from meters
    % (we estimate from 45 deg azimuth)
    if num_p > 1
        [lon2,lat2] = m_fdist_par(lon_g,lat_g,45,hh_m,num_p);
    else
        [lon2,lat2] = m_fdist(lon_g,lat_g,45,hh_m);
    end
    if ~isempty(find(isnan(lon2), 1))
       lon2 = fillmissing(lon2,'linear');
       lat2 = fillmissing(lat2,'linear');
    end
    % switch back to -180 to 180
    lon2(lon2 > 180) = lon2(lon2 > 180) - 360;
    % get the actual distance
    hh_m = sqrt((lon2 - lon_g).^2 + (lat2 - lat_g).^2);
    
    % Get min of all the criteria
    if dist_param > 0 && ~ocean_only
        hh_m = min(hh_m,hh(:,:,1));
    end
else
    hh_m = squeeze(hh);
end

% min edgelength
hh_m(hh_m < min_el) = min_el;
% max edgelength
hh_m(hh_m > max_el) = max_el;

% % relax gradient with Mesh2D's hill-climb algorithm...
% % modified to elminate the need for an unstructured mesh...
%hfun = reshape(hh_m,[],1);
hfun = zeros(size(hh_m,1)*size(hh_m,2),1); 
nn = 0;
% renumber so row varies the fastest 
for ipos = 1 : size(lon_g,1)
    for jpos = 1 : size(lon_g,2)
        nn = nn + 1;
        hfun(nn,1) = hh_m(ipos,jpos);
    end
end
eglen = abs(lon_g(2)-lon_g(1)); dhdx = dist_param; 
ny = size(lon_g,2);
disp('Relaxing the gradient...');
tic
[hfun,flag] = limgradStruct(ny,eglen,hfun,dhdx,sqrt(length(hfun)));
toc
if flag == 1
    disp('Gradient relaxing converged!');
else
    disp(['Gradient relaxing did not converge, '
          'please check your edge functions (run with plot_on == 2 or 3)']);
end
%hh_m = reshape(hfun,[],ny);
% % reshape it back 
nn = 0;
for ipos = 1 : size(lon_g,1)
    for jpos = 1 : size(lon_g,2)
        nn = nn+1;
        hh_m(ipos,jpos) = hfun(nn);
    end
end

% Plotting the overall edgelength function
if plot_on >=1
    figure;
    m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],...
                      'lat',[bbox(2,1) bbox(2,2)]);
    m_pcolor(lon_g,lat_g,hh_m); shading interp;
    m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
    title('Press 1 to accept, 0 to abort...'); colorbar;
    ierr = abort; if ierr == 1; return; end
    print('Overall','-dpng','-r300')
end

% Make the overall interpolant
disp('Building the finalized edge function interpolant...')
F = griddedInterpolant(lon_g,lat_g,hh_m,'linear');
clear hh hh_m hfun2 hfun3 bathy

%% Call the distmesh iteration routine
fprintf(1,' ------------------------------------------------------->\n') ;
disp('Entering DistMesh2d');
close all;
[p,t] = distmesh2d(@fd,@fh,min_el,bbox',ini_p,fix_p,itmax,plot_on,nscreen);
disp('Exiting Distmesh2d.m');
fprintf(1,' ------------------------------------------------------->\n') ;

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
mq_l = min(tq.qm);
disp(['number of nodes is ' num2str(length(p))])
disp(['mean quality is ' num2str(mq_m)])
disp(['min quality is ' num2str(mq_l)])  
%
if plot_on >= 1
    simpplot(p,t); title('Finalized');
end
fprintf(1,' ------------------------------------------------------->\n') ;
disp('Exiting Distmesh_prep...');
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
    end

    function d = fd( p, instance )
        % The distance function
        if instance == 0
            % used for kicking points out of the box.
            % mdl0 contains the ocean boundary
            if ~isempty(meshfile)
                d = dpoly_fp2(p,mdl0,polygon,floodplain_polygon,Fb,[]);
            else
                poly = polygon_struct.outer;
                if ~isempty(polygon_struct.inner)
                    poly = [poly; NaN NaN; polygon_struct.inner];
                end
                if isempty(dOB)
                    d = dpoly(p,poly,mdl0,num_p);
                else % reduce number of pts checked by only checking pts nearby (dclose) outer boundary
                    dclose = 2*dxx; % the distance from boundary that we check points
                    np = length(p);
                    ptest(:,1:2) = p; ptest(:,3) = 1:length(p);
                    ptest(FOB(ptest(:,1),ptest(:,2)) < -dclose,:) = [];
                    dtest = dpoly(ptest(:,1:2),poly,mdl0,num_p);
                    d = ones(np,1)*-999; % all points are assumed in
                    d(ptest(:,3)) = dtest; % except for dtest
                    clearvars ptest psave dtest
                end
            end
        elseif instance == 1
            % used for the creation of the distance function (sans ocean bou).
            % mdl1 does not contain the ocean boundary
            if ~isempty(meshfile)
                d = dpoly_fp2(p,mdl1,polygon,floodplain_polygon,[],[]);
            else
                poly = polygon_struct.outer;
                if ~isempty(polygon_struct.inner)
                    poly = [poly; NaN NaN; polygon_struct.inner];
                end
                d = dpoly(p,poly,mdl1,num_p);
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

    function ierr = abort( )
        go_on  = 0;
        go_on = input('Press 1 to go on, 0 to abort...');
        if go_on ~= 1
            p = []; t = [];
            ierr = 1; 
            return
        else
            ierr = 0; 
            return;
        end
    end

    function dx = m_idist_par(lon1,lat1,lon2,lat2,num_p)
        if isempty(gcp)
            parpool('local',num_p);
        end
        Pool = gcp('nocreate');
        dx  = zeros(size(lon1));
        np1 = size(lon1,1);
        for idx = 1:num_p
            ns1 = int64((idx-1)*np1/num_p)+1;
            ne1 = int64(idx*np1/num_p);
            f(idx) = parfeval(Pool,@m_idist,1,...
                              lon1(ns1:ne1,:),lat1(ns1:ne1,:),...
                              lon2(ns1:ne1,:),lat2(ns1:ne1,:));
        end
        for idx = 1:num_p
            [idx_t, dx_t] = fetchNext(f); % Get results into a cell array
            ns1 = int64((idx_t-1)*np1/num_p)+1;
            ne1 = int64(idx_t*np1/num_p);
            dx(ns1:ne1,:) = dx_t;
        end
    end

    function [lon2,lat2] = m_fdist_par(lon1,lat1,angle,hh_m,num_p)
        if isempty(gcp)
            parpool('local',num_p);
        end
        Pool = gcp('nocreate');
        np1  = size(hh_m,1) ;
        lat2 = zeros(size(hh_m));
        lon2 = zeros(size(hh_m));
        for idx = 1:num_p
            ns1 = int64((idx-1)*np1/num_p)+1;
            ne1 = int64(idx*np1/num_p);
            f(idx) = parfeval(Pool,@m_fdist,2,...
                              lon1(ns1:ne1,:),lat1(ns1:ne1,:),...
                              angle,hh_m(ns1:ne1,:));
        end
        for idx = 1:num_p
            [idx_t, lon2_t, lat2_t] = fetchNext(f); % Get results into a cell array
             ns1 = int64((idx_t-1)*np1/num_p)+1;
             ne1 = int64(idx_t*np1/num_p);
             lon2(ns1:ne1,:) = lon2_t;
             lat2(ns1:ne1,:) = lat2_t;
        end
    end

    function [outer,ierr] = cornerify(outer,boubox,plot_on)
        ierr = 0;
        if plot_on < 1
           fprintf(1,' -------------------------------------------------------\n') ;...
           disp('FATAL: You have provided a segment, you must first run in interactive mode...');
           p = []; t = [];....
           fprintf(1,' -------------------------------------------------------\n') ;
           ierr = 1;
           return;
        end
        fprintf(1, [ ...
            ' INFO: You have provided a segment as your contourfile. We must \n', ...
            ' add two corners to turn it into a closed polygon \n']);
        
        title('Please close the outer segment to form a polygon...')
        % Here, we find the portion of the outer polygon that is outside the
        % bbox. Then we find the maximum and minimum coordinates of the
        % outer segment not inside the bbox and extend them to use as corners to form a
        % rectangle that closes our outer segment.
        outer(end,:) = []; % remove trailing NaN
        kept = outer(inpoly(outer,boubox)~=1,:); % find outer segment outside of bbox
        m_plot(boubox(:,1),boubox(:,2),'r-'); hold on;
        m_hatch(boubox(:,1),boubox(:,2),'single',30,5,'color','k'); % ...with hatching added.
        
        xlen = max(diff(boubox),1);
        side = input('PROMPT: Which side (right=1,left=0) of the coastline is the ocean on?');
        
        botcorner = min(kept); botcorner(1) = botcorner(1) + xlen*(-1)^side;
        m_plot(botcorner(1),botcorner(2),'ro');
        
        topcorner =  [botcorner(1),botcorner(2) + ...
                                   abs(max(outer(:,2)) - min(outer(:,2)))];
        m_plot(topcorner(1),topcorner(2),'ro');
        
        outer = [outer; botcorner];
        outer = [outer; topcorner];
        outer = [outer; max(kept)];
        outer = [outer; NaN NaN];
    end
end
