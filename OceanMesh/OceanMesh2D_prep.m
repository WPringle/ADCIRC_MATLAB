function [p,t] = OceanMesh2D_prep(contourfile,bathyfile,feature_vecs,...
     bbox,min_el,max_el,max_el_ns,dist_param,wl_param,...
    slope_param,channel_param,dt,bounds,itmax,...
    plot_on,ini_p,fix_p,nscreen,num_p,edgefx,DriverCounter,outfiname)
% Builds the edge function and necessary dependencies for the calling mesh
% smoothing algorithm.
%  V0 : Initial development by William Pringle, Damrongsak Wirasaet, 2016-Feb 2017
% V1:  Combining floodplain+coastal meshing,overall improvements+bugfixes, and
%      polygonal selection tool for floodplain meshing by KJR 2017-April-2017-June
% V2.0 C++ ANN KD-Tree for NN-searches, feature size edgefunction, and general improvements by KJR July 2017.
% V2.x Making bathyfiles, feature size and CFL limiter options more general, small changes in distmesh2d. WJP Aug 2017
% V2.5 Can loop through more than one dataset,improvements to DEM handling and channel scale edge function, improvements to the feature size function.
%      improvements to clipping algorithm to form outer-polygon (now using SutherlandHodgman algorithm).
%      Corrections in the fd (distance function) to avoid problems with inpoly. KJR Aug 2017.
% V3.0: Combining coastal + floodplain functionality with the idea of feature vectors. Deprecating uncessary options. KJR OCTOBER 2017.
%
%
% INPUTS:  
%          contourfile: filename(s) of shape file(s) (must be a cell)
%          bathyfile  : filename(s) of the bathymetric dataset in NETCDF
%                       format
%         feature_vecs: A cell-array of shapefiles containing features that you
%                       want to preserve in the triangulation.
%          bbox       : the bounding box to mesh within in following format:
%                       [minlon maxlon;
%                        minlat maxlat];
%          min_el     : The minimum (and initial) edgelength of the mesh
%          max_el     : The maximum edgelength of the mesh.
%          max_el_ns  : The maximum edgelength of the mesh on the coastline.
%          dist_param : The fraction that the edgelength should change
%                       with distance. Set 0 to ignore
%          wl_param   : Parameter in wavelength function, (set 0 to ignore)
%                       edgelength = T/wl_param*sqrt(gH), where T is the
%                       period of the M2 wave and H is the bathymetric depth
%          slope_param: Parameter in slope function, (set 0 to ignore)
%                       edgelength = 2*pi/slope_param*H/abs(grad(H)),
%                       where H is the bathymetric depth
%          channel_parm: Parameter that determines how much resolution is
%                        used in areas of strong gravity driven flow (normally 0.5e-6)
%                        or large upslope area
%          dt         : Timestep to use for CFL limiter. dt < 0 to ignore
%                       If dt = 0 , automatic based on distance function
%          bounds     : Where meshfile is present (to mesh floodplain)
%                       bounds(1) is the upper bound of the topo to mesh
%                       within and bounds(2) is the max dist from coastline
%                       to mesh within
%          itmax      : Maximum number of iterations in distmesh function
%          plot_on    : choose whether to plot map and mesh during calc or
%                       not (Yes = 1, No = 0)
%          ini_p      : optional initial distribution of p
%          fix_p      : optional fixed node points
%          nscreen    : Frequency of iterations to save a temporary mesh
%                       and display the mesh quality
%          num_p      : number of parallel processors required(<=1 for serial)
%          edgefx     : Edge function interpolant. A gridded interpolant from a previous 
%	               		go at the same problem.
%
% OUTPUTS: p          : the node positions in lon,lat
%          t          : the triangulation of the mesh (in counter-clockwise order).
%
%
%
% Reference for wavelength and slope function:
% Lyard, F., Lefevre, F., Letellier, T., & Francis, O. (2006).
% Modelling the global ocean tides: modern insights from FES2004.
% Ocean Dynamics, 56(5�6), 394�415. http://doi.org/10.1007/s10236-006-0086-x
%
% Reference for feature size function:
% Koko J., A Matlab mesh generator for the two-dimensional finite element method, 
% Applied Mathematics and Computation 250, p. 650-664 (215)
% 
% Reference for Topotoolbox
% Schwanghart, W., Scherler, D. (2014): TopoToolbox 2 – MATLAB-based software 
% for topographic analysis and modeling in Earth surface sciences. Earth Surface Dynamics, 2, 1-7.
% [DOI: 10.5194/esurf-2-1-2014, Discussion paper: 10.5194/esurfd-1-261-2013]. 
% open access and software available here [DOI: 10.1594/IEDA/100175].
%
%% Build DEM
% Read in the GeoTIFF/ERSI ASCII text file as a DEM object
disp('INFO: Building DEM object...');
if exist(bathyfile{DriverCounter}, 'file') == 2
    DEM           = GRIDobj(bathyfile{DriverCounter});
    DEMc          = DEM.crop(bbox(1,:),bbox(2,:));
    clearvars DEM                                                          % release full DEM since no longer necessary after cropping.
    info(DEMc)                                                             
    [lon,lat]     = DEMc.getcoordinates;
    lat           = flipud(lat);
    [lon_g,lat_g] = ndgrid(lon(:),lat(:));                                 % must be in ndgrid format for gridded interpolant
    bathy         = DEMc.Z;
    bathy         = flipud(bathy)';
    bathy(bathy==-9999)=NaN;
    bathyres      = DEMc.cellsize;
    Fb = griddedInterpolant(lon_g,lat_g,bathy,'linear','nearest'); 
else
    disp('INFO: Bathyfile does not exist...');
end

% fill in NaNs with srtm15+,
if ~isempty(find(isnan(bathy), 1))
    if length(bathyfile) > 1
        disp('INFO: Filling NaNs from bathyfile fix...')
        lond = double(ncread(bathyfile{end},'lon'));
        latd = double(ncread(bathyfile{end},'lat'));
        Id = find(lond > bbox(1,1) & lond < bbox(1,2));
        Jd = find(latd > bbox(2,1) & latd < bbox(2,2));
        lond = lond(Id); latd = latd(Jd);
        % read only part of the DEM necessary
        ncid = netcdf.open(bathyfile{end},'NC_NOWRITE');
        bathyFIX = double(netcdf.getVar(ncid,2,[Id(1) Jd(1)],[length(Id) length(Jd)]));
        netcdf.close(ncid)
        [lon_gd,lat_gd] = ndgrid(lond,latd);
        Fb_fix = griddedInterpolant(lon_gd,lat_gd,bathyFIX);
        bathy(isnan(bathy)) = Fb_fix(lon_g(isnan(bathy)),lat_g(isnan(bathy)));
        clearvars lond latd Id Jd bathyFIX lon_gd lat_gd Fb_fix
    else
        disp('WARNING: Some NaNs in the bathy but no secondary bathy provided')
    end
end
%bathy = -bathy;                                                            %--> is the bathy negative below the geoid or positive? 
%% Process shapefile
if ~isempty(contourfile) 
    disp('Operating in coastal mode...');
    % Make the bounding box 5 x 2 matrix in clockwise order
    boubox = [bbox(1,1) bbox(2,1); bbox(1,1) bbox(2,2); ...
        bbox(1,2) bbox(2,2); bbox(1,2) bbox(2,1); ...
        bbox(1,1) bbox(2,1)];
    polygon_struct = [];
    polygon        = [];
    
    disp('Reading in shapefile...')
    % Read polygon from shape file, and apply spacing interpolant
    polygon_struct = Read_shapefile( contourfile, polygon, bbox, ...
        abs(min_el), 0 );
    
    % Smooth the polygons and segments that exist using a 5-point moving average
    poly_count = 0; oc = 0; inr = 0; ml = 0; ocean_only = 0;
    if ~isempty(polygon_struct.outer)
        oc = oc + 1; poly_count = poly_count + 1;
        [polygon_struct.outer,closed] = smooth_coastline(polygon_struct.outer,5,0);
    end
    if ~isempty(polygon_struct.mainland)
        ml = ml + 1; poly_count = poly_count + 1;
        polygon_struct.mainland = smooth_coastline(polygon_struct.mainland,5,0);
    end
    if ~isempty(polygon_struct.inner)
        inr = inr + 1; poly_count = poly_count + 1;
        polygon_struct.inner = smooth_coastline(polygon_struct.inner,5,0);
    end
    % handle the case when no segments are present are in the bbox.
    if oc == 1 && poly_count == 1
        ocean_only = 1;
        disp('INFO: Meshing the open ocean, no land or island segments present within bbox...');
    end
    
    if plot_on >= 1
        bufx = 0.2*(bbox(1,2) - bbox(1,1));
        bufy = 0.2*(bbox(2,2) - bbox(2,1));
        m_proj('Mercator','long',[bbox(1,1) - bufx, bbox(1,2) + bufx],...
            'lat',[bbox(2,1) - bufy, bbox(2,2) + bufy])
        if(ml~=0)
            m_plot(polygon_struct.mainland(:,1),polygon_struct.mainland(:,2),...
                'm-','linewi',2); hold on;
        else
            m_plot(polygon_struct.inner(:,1),polygon_struct.inner(:,2),...
                'm-','linewi',2); hold on;
        end
        m_plot(boubox(:,1),boubox(:,2),'k-o','linewi',2,'MarkerSize',15);
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    end
    
    % Cornerify: create artifical outer bounding box to turn segment into
    % polygon. Lets test for this automatically
    if ~closed && ~ocean_only
        [temp] = cornerify(polygon_struct.mainland,boubox);
        [la,lo]=interpm(temp(:,2),temp(:,1),abs(min_el)/2);
        polygon_struct.outer = [];
        polygon_struct.outer(:,1) = lo; polygon_struct.outer(:,2) = la;
    end
    if plot_on >=1
        close all;
        m_proj('Mercator','long',[bbox(1,1) - bufx, bbox(1,2) + bufx],...
            'lat', [bbox(2,1) - bufy, bbox(2,2) + bufy])
        hold on;
        % Plot the polygons that exist
        if inr == 1
            m_plot(polygon_struct.inner(:,1),polygon_struct.inner(:,2),'m');
        end
        if ml == 1
            m_plot(polygon_struct.mainland(:,1),polygon_struct.mainland(:,2),'k-')
        end
        m_plot(boubox(:,1),boubox(:,2),'k','linewi',2);
        m_hatch(polygon_struct.outer(:,1),polygon_struct.outer(:,2),'single',30,5,'color','k'); % ...with hatching added.
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        title('Meshing the union of the filled areas...press 1 to accept or 0 to abort...');
        save(['PG_',outfiname,'.mat'], 'polygon_struct','ocean_only' ,'inr' ,'oc','ml', 'poly_count')
        if(~isempty(fix_p))
            save(['FixP_',outfiname,'.mat'],'fix_p');
        end
        ierr = abort(1); if ierr==1; return; end
    end
elseif ~isempty(feature_vecs)
    disp('Operating in feature vector mode...');
    % If we don't have an outer polygon already then make it by bbox
    polygon_struct = []; polygon=[];
    polygon_struct.outer = [bbox(1,1) bbox(2,1);
        bbox(1,1) bbox(2,2);
        bbox(1,2) bbox(2,2);
        bbox(1,2) bbox(2,1);
        bbox(1,1) bbox(2,1)];
    [la,lo]=interpm(polygon_struct.outer(:,2),polygon_struct.outer(:,1),abs(min_el)/2);
    polygon_struct.outer=[]; 
    polygon_struct.outer(:,1)=lo;
    polygon_struct.outer(:,2)=la; 
    polygon_struct.outer=[polygon_struct.outer;NaN,NaN]; 
    % the first feature_vec is the 0-m contour, make this the mainland 
    dmy  = Read_shapefile( {feature_vecs{1}}, polygon, bbox, ...
        abs(min_el), 0 );
    polygon_struct.mainland = dmy.mainland;
    polygon_struct.inner    = dmy.inner;
    ocean_only = 0; inr = 0; oc = 0; ml = 0; poly_count = 0; 
    save(['PG_',outfiname,'.mat'], 'polygon_struct','ocean_only' ,'inr' ,'oc','ml', 'poly_count')
    if(~isempty(fix_p))
        save(['FixP_',outfiname,'.mat'],'fix_p');
    end
else
    % you have already run in interactive mode just load the polygon_struct
    disp(['INFO: Neither the contourfile or the meshfile were provided, '...
        'trying to load in the polygon_struct from a previous interactive session...']);
    if exist(['PG_',outfiname,'.mat'],'file')
        disp(['INFO: Loading in polygon_structure: PG_',outfiname,'.mat'])
        load(['PG_',outfiname,'.mat'])
    else
        disp('FATAL: polygon_struct does not exist, exiting...')
        return;
    end
    if exist(['FixP_',outfiname,'.mat'],'file')
        disp(['INFO: Loading in Fixed points: FixP_',outfiname,'.mat'])
        load(['FixP_',outfiname,'.mat'])
    end
end
%% Build KD-Trees
% For the distance function used to calculate edgelength function
if ~ocean_only
    disp('Building KD-Tree (mdl1) with mainland and island segments...');
    c_pts = [polygon_struct.inner; polygon_struct.mainland];
    c_pts(isnan(c_pts(:,1)),:) = [];
    mdl1  = ann(c_pts');
else
    mdl1 = [];
end

% For distance function used to determine whether points and elements
% are inside the region you want to mesh
disp('Building KD-Tree (mdl0) with ocean boundary, mainland, and island segments...');
if ~ocean_only && isempty(feature_vecs)
    poly = [polygon_struct.outer; NaN NaN; ...
        polygon_struct.inner];
    poly(isnan(poly(:,1)),:) = [];
    mdl0 = ann(poly');
    poly = [];
else
    poly = polygon_struct.outer;
    poly(isnan(poly(:,1)),:) = [];
    mdl0 = ann(poly');
    poly = [];
end
%% Build Edge Fx
if ~isempty(edgefx)    
    % The edge function shall be named "hh_m" and have nx by ny by 4
    % dimensions and be in degrees.
    disp('Loading edge function...');
    tic
    load(edgefx);
    toc
else
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
        d = reshape(d,size(lon_g,1),[]);
        nearshore = d > -0.01 & d < eps;
        if min_el > 0
            disp('   using feature size function...');
            % calculate the gradient of the distance function.
            [ddx,ddy] = gradient(d,bathyres);
            d_fs = sqrt(ddx.^2 + ddy.^2);
            % find singularties in the distance function that are
            % within the poly to get the medial axis.
            % Lets only create medial points in places that are
            % sufficiently far away from the coastline to capture the
            % feature and thus aren't suprious.
            d_fs = reshape(d_fs,[],1); d = reshape(d,[],1);
            x_kp = x_v(d_fs < 0.90 & d < -eps); %-2*bathyres);
            y_kp = y_v(d_fs < 0.90 & d < -eps); %-2*bathyres);
            if(plot_on >=2)
                hold on; m_plot(x_kp,y_kp,'r.');% plot the medial points
            end
            % Now get the feature size along the coastline
            [~, dPOS] = WrapperForKsearch([x_kp,y_kp]',[x_v,y_v]');
            % reshape back
            d = reshape(d,size(lon_g,1),[]); dPOS = reshape(dPOS,size(lon_g,1),[]);
            % feature_size is distance from medial axis plus distance to
            % coastline. We put the dist_param on d to make size larger
            % from coastline. min_el is then feature_size*2/R where R is
            % number of elements to model the feature
            fs = (2*(dPOS-d))/3;
            fs(nearshore & fs >  max_el_ns)=  max_el_ns; %--enforce near shore max ele
            hh(:,:,nn) = fs; %(2*(dPOS-d))/3;
            clear x_kp y_kp d_fs dPOS
        elseif min_el < 0
            disp('   not using feature size function...');
            % add into edge function
            hh(:,:,nn) = abs(min_el) - dist_param*d ;
        end
        
        % Plot the edgelength function for distance
        if plot_on >= 2
            figure;
            m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)])
            m_pcolor(lon_g,lat_g,hh(:,:,nn)); shading interp
            hold on; m_plot(polygon_struct.outer(:,1),polygon_struct.outer(:,2),'k-','linewi',1);
            cb = colorbar; ylabel(cb,'edgelength in degrees');
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
            title('Distance function');
            ierr = abort(1); if ierr == 1; return; end
            print('Distancefunction','-dpng','-r300')
        end
        
    end
    
    % Wavelength...
    if wl_param > 0
        disp('   Building the wavelength edge function...');
        nn = nn + 1;
        g = 9.807;
        period = 12.42*3600; % M2 period in seconds
        temp       = period*sqrt(g*abs(bathy))/wl_param;
        temp(bathy > eps) = NaN; %--> above sea level is not defined.
        hh(:,:,nn) = temp;
       
        if plot_on >= 2
            figure;
            m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)]);
            m_pcolor(lon_g,lat_g,real(hh(:,:,nn))/1000); shading interp
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
            title('Wavelength of the M2 edge function'); cb=colorbar;ylabel(cb,'km');
            ierr = abort(1); if ierr == 1; return; end
            print('Wavelengthfunction','-dpng','-r300')
        end
    end
    
    % Slope...
    if slope_param > 0
        disp('   Building the slope edge function...');
        nn        = nn + 1;
        % Converting lon and lat to metres so we can evaluate the slope
        lon2 = single([lon_g(2:end,:); lon_g(end,:)+(lon(end)-lon(end-1))]);
        lat2 = single([lat_g(:,2:end) lat_g(:,end)+(lat(end)-lat(end-1))]);
        % lon
        if num_p > 1
            dx = m_idist_par(lon_g,lat_g,lon2,lat_g,num_p*8);
        else
            dx = m_idist(lon_g,lat_g,lon2,lat_g);
        end
        % lat
        if num_p > 1
            dy = m_idist_par(lon_g,lat_g,lon_g,lat2,num_p*8);
        else
            dy = m_idist(lon_g,lat_g,lon_g,lat2);
        end
        % Get the coordinate matrices for dx and dy
        dx = cumsum(dx,1); dx = cumsum(dx,2);
        dy = cumsum(dy,1); dy = cumsum(dy,2);
        
        % Calculate the gradient of bathy
        [b_y,b_x] = gradient(bathy,dy,dx); clearvars dx dy
        b_slope = sqrt(b_x.^2 + b_y.^2); clearvars b_x b_y
        % slope func
        hh(:,:,nn) = abs(2*pi*bathy./b_slope/slope_param);
        clear lon2 lat2 b_slope
        clearvars dx dy d nv k
    end
    
    % channel scale
    if(channel_param > 0)
        disp('   Building the confluence scale edge function...');
        nn=nn+1;
        DEMc          = DEMc.elevateminima;                                 % get rid of regional sinks
        FD            = FLOWobj(DEMc,'preprocess','carve','mex',true);      % create flow obj and carve dem to fill in local sinks.
        A             = flowacc(FD);                                        % n cells draining into cell ie SFD algo.
        A             = dilate(A,ones(3));                                  % enlarges the features to better capture on the scales we're interested in.
        channels = abs(DEMc.Z)/A;                                           % volume of cell divided by upslope area
        channel_norm = (min_el*111000)/channel_param ;                      % numerator is the mininum element size in meters. Demoninator has to be representative of flow scale that gets min_el.
        hh(:,:,nn)    = flipud(channels.Z)'*channel_norm;                   % normalize it so the min_el corresponds to a specific scale.
        
        if plot_on >= 2
            figure;
            m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],'lat',[bbox(2,1) bbox(2,2)]);
            m_pcolor(lon_g,lat_g,real(hh(:,:,nn))); shading interp
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
            title('Confluence scale edge function'); cb=colorbar;ylabel(cb,'km');
            ierr = abort(1); if ierr == 1; return; end
            print('ChannelFunction','-dpng','-r300')
        end
        clearvars lon_g2 lat_g2 channel_scale_grd FD A
    end
    
    %% EnforceBounds
    % Get min of slope, wavelength and channel scale, setting equal to hh_m
    if wl_param > 0 || slope_param > 0 || channel_param > 0
        if dist_param > 0 && ~ocean_only % dist param is already in degrees
            hh_m = min(hh(:,:,2:end),[],3);
            hh(:,:,2:end)= []; % these are large arrays-release the memory
        else
            hh_m = min(hh(:,:,1:end),[],3);
        end
        dlon=cosd(lat)'; % distance between longitude's
        dlon=repmat(dlon,[size(lon_g,1),1]);
        hh_m=(abs(dlon)./111000).*hh_m;
        
        % Get min of all the criteria
        if dist_param > 0 && ~ocean_only
            hh_m = min(hh_m,hh(:,:,1));
        end
    else
        hh_m = squeeze(hh);
    end
    clear lat2 lon2 x_v y_v% these can be very large arrays, need to release
    
    % enforce min edgelength
    hh_m(hh_m < abs(min_el)) = abs(min_el);
    % enforce max edgelength
    hh_m(hh_m > max_el) = max_el;
    
    % relax gradient with Mesh2D's hill-climb algorithm...
    % modified to elminate the need for an unstructured mesh...
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
    [hfun,flag] = limgradStruct(ny,eglen,hfun,dhdx,sqrt(length(hfun)));
    if flag == 1
        disp('Gradient relaxing converged!');
    else
        disp(['Gradient relaxing did not converge, '
            'please check your edge functions (run with plot_on == 2)']);
    end
    % reshape it back
    nn = 0;
    for ipos = 1 : size(lon_g,1)
        for jpos = 1 : size(lon_g,2)
            nn = nn+1;
            hh_m(ipos,jpos) = hfun(nn);
        end
    end
    
    % Limit CFL if dt >= 0, dt = 0 finds dt automatically
    if dt >= 0
        g = 9.807; descfl = 0.25;
        if dt == 0
            % Find min allowable dt based on distance function
            if dist_param > 0 && ~ocean_only
                hh_d = hh(:,:,1); hh_d(hh_d <= 0) = NaN;
                hh_d(hh_d < min_el ) = min_el;
                [~,loc] = max(hh_d(:));
                % This ensures only positive depths taken into account
                bb = max(sign(bathy(loc))*bathy,0);
                dt = descfl*hh_d*111000./sqrt(g*bb);
                [dt, ~] = min(dt(:));
                clear hh_d
            else
                disp('Error: cannot use automatic CFL limiter with no dist_param or no land boundaries')
                abort(0);
            end
        end
        disp(['Enforcing CFL condition of ',num2str(descfl),...
            ' for an expected simulation dt of ',num2str(dt)]);
        cfl = dt*sqrt(g*abs(bathy))./(hh_m*111000); %--this is your cfl
        dxn = sqrt(g*abs(bathy))*dt/descfl;      %--assume simulation time step of dt sec and cfl of dcfl;
        hh_m( cfl > descfl ) = dxn( cfl > descfl )/111000;   %--in degrees
        clear cfl dxn
    end
    clear bathy hh
    
    % Plotting the overall edgelength function
    if plot_on >=1
        figure;
        m_proj('Mercator','long',[bbox(1,1) bbox(1,2)],...
            'lat',[bbox(2,1) bbox(2,2)]);
        m_pcolor(lon_g(1:10:end,1:10:end),lat_g(1:10:end,1:10:end),hh_m(1:10:end,1:10:end)); shading interp;
        m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
        title('Press 1 to accept, 0 to abort...'); colorbar;
        ierr = abort(1); if ierr == 1; return; end
        print('Overall','-dpng','-r300')
    end
    
    % Make the overall interpolant
    disp('Building the finalized edge function interpolant...')
    F = griddedInterpolant(lon_g,lat_g,hh_m,'linear','nearest');
    save -v7.3 EdFxnIntrp F
    clear hh_m hfun2 hfun3
end %--end load edgefx

% If present, load feature vectors and store as fixed points.
if ~isempty(feature_vecs)
    disp('Reading in feature vectors...');
    count = 0;
    for fname=feature_vecs
        % Read the structure
        S = shaperead(fname{1},'BoundingBox',bbox);
        
        % Get rid of unwanted components
        Fi = fieldnames(S);
        D = struct2cell(S);
        S = cell2struct(D(3:4,:), Fi(3:4));
        
        if ~isempty(S)
            % Keep the following polygons
            for i = 1 : length(S)
                for j = 1 : length(S(i).X)-1
                    count = count + 1;
                    fix_p(count,:) = [S(i).X(j),S(i).Y(j)];
                end
            end
        end
    end

end

%% MeshSmoothing
% Call the modified distmesh smoothing routine
fprintf(1,' ------------------------------------------------------->\n') ;
disp('Entering DistMesh2d');
close all;
[p,t] = distmesh2d_plus(@fd,@fh,abs(min_el),bbox',ini_p,fix_p,itmax,plot_on,nscreen);
disp('Exiting Distmesh2d.m');
fprintf(1,' ------------------------------------------------------->\n') ;


%% PostProcessMesh
%% Remove portion of the mesh that exceeds your distance AND depth criteria
if ~isempty(feature_vecs)
    disp('Trimmming mesh...')
    [p,t] = trim_mesh(p,t);
    [p,t] = Fix_bad_edges_and_mesh(p,t,0); 
end

% Fixing up the mesh automatically
if isempty(fix_p)
    disp('Check for and fix poor quality triangles...');
    [p,t] = Fix_bad_edges_and_mesh(p,t,0);
    tq=gettrimeshquan(p,t); t(abs(tq.qm)<0.10,:) = [];
    [p,t] = direct_smoother_lur(p,t);
    [p,t] = Fix_bad_edges_and_mesh(p,t,0);
end
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
        pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;      % Compute centroids
        dis = fd(pmid,2);
        t(dis > bounds(2) && Fb(pmid) > bounds(1),:)=[];   % Remove triangles > than maxDist away from coast AND greater than MaxEle above MSL.
        [p,t] = fixmesh(p,t);                              % Clean up disjoint vertices.
    end


    function d = fd( p, instance )
        % The distance function
        if instance == 0
            % used for kicking points out of the box.
            % mdl0 CONTAINS the ocean boundary
            poly = polygon_struct.outer;
            d = dpoly(p,poly,mdl0,num_p,poly);
            % IF OUTSIDE BUT APPEARS INSIDE
            bad = find((p(:,1) < bbox(1,1) | p(:,1) > bbox(1,2) | ...
                p(:,2) < bbox(2,1) | p(:,2) > bbox(2,2)) & d < 0);
            
            if ~isempty(bad)
                d(bad) = -d(bad);
            end
            
        elseif instance == 1
            % used for the creation of the distance function (sans ocean bou).
            % mdl1 does NOT contain the ocean boundary
            poly = polygon_struct.outer;
            if ~isempty(polygon_struct.inner)
                poly = [poly; NaN NaN; polygon_struct.inner];
            end
            d = dpoly(p,poly,mdl1,num_p,[polygon_struct.inner; polygon_struct.mainland]);
            % IF INSIDE BUT APPEARS OUTSIDE
            bad = find((p(:,1) >= bbox(1,1) | p(:,1) <= bbox(1,2) | ...
                p(:,2) >= bbox(2,1) | p(:,2) <= bbox(2,2)) & d > 0);
            
            if nnz(bad) > 0
                d(bad) = -d(bad);
            end
        end
        return ;
    end


    function h = fh( p )
        % The edge length function, desides on the edgelength at this
        % location, based on distance, wavelength, slope, and feature size
        % fxs
        h = F(p);
        return ;
    end

    function ierr = abort(go_on)
        go_on  = 1;
        go_on = input('Press 1 to go on, 0 to abort...'); %- comment this out if you don't want the program to wait for you!
        if ~go_on
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

    function [latout,lonout]=densify(latin,lonin)
        % Densify the outer polygon (fills gaps larger than half min edgelength).
        % TODO REWRITE THIS SO THE USER DOES NOT NEED LICESNSE 
        [latout,lonout] = interpm(latin,lonin,abs(min_el));
        return;
    end

    function [clippedPolygon] = cornerify(segment,boubox)
        % INPUTS: segment: A segment/polygon or many segments or many polygons stored as an n x 2 array with the x,y coordinates--last entry are NaNs
        %         boubox:  An 5x2 array containing the x,y locations of the
        %         corners of the region you want to mesh.
        % OUTPUTS: An array of n x 2 size containing the x,y coordinates of the
        %          region that represents the extent of the intersection between segment and boubox.
        
        %--INFO: This algorithm closes segment(s) to form polygon(s) so we can
        % perform polygon clipping using the Sutherland-Hodgman algorithm.
        % This gives an outer polygon for the distmesh algorithm to work.
        % Sutherland-Hodgman algorithm source:
        % https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#MATLAB_.2F_Octave
        %--kjr 2017, UND.
        %
        lidx   = find(isnan(segment(:,1))); %--find segments.
        lidx   = [1;lidx];
        outer  = []; clippedPolygon = [];
        disp('There are some segments that must be closed. We will display the segments on a map, then follow the prompt...');
        
        for noSegs = 1 : length(lidx)-1
            close all;
            if(noSegs==1)
                st = 1;
            else
                st = lidx(noSegs)+1;
            end
            ed = lidx(noSegs+1);
            % if a poly, then add it to outer but don't try to close it.
            if segment(st,1) == segment(ed-1,1) || segment(st,2) == segment(ed-1,2)
                clippedPolygon = [clippedPolygon; segment(st:ed-1,:); NaN,NaN];
                disp('Found a polygon, continuing...');
                fprintf(1,' ------------------------------------------------------->\n') ;
                fprintf(1,' ------------------------------------------------------->\n') ;
                continue
            end
            tbbox = [min(segment);max(segment)]';
            tbufx = 0.2*(tbbox(1,2) - tbbox(1,1));
            tbufy = 0.2*(tbbox(2,2) - tbbox(2,1));
            m_proj('Mercator','long',[tbbox(1,1) - tbufx, tbbox(1,2) + tbufx],...
                'lat',[tbbox(2,1) - tbufy, tbbox(2,2) + bufy]);
            m_plot(boubox(:,1),boubox(:,2),'k','linewi',2);
            m_grid('xtick',10,'tickdir','out','yaxislocation','left','fontsize',7);
            

            hold on; m_text(segment(st,1),segment(st,2),'Start','FontSize',15);
            hold on;m_plot(segment(st,1),segment(st,2),'rx','MarkerSize',15);
            hold on;m_text(segment(ed-1,1),segment(ed-1,2),'End','FontSize',15);
            hold on;m_plot(segment(ed-1,1),segment(ed-1,2),'rx','MarkerSize',15);
            hold on;m_plot(segment(:,1),segment(:,2),'k-');
            hold on;m_plot(segment(st:ed-1,1),segment(st:ed-1,2),'b','linewi',2);
            for corn = 1 : 4
                hold on;m_text(boubox(corn,1)-0.001,boubox(corn,2)+0.001,num2str(corn),'fontsize',15);
            end
            fprintf(1,[...
                ' INFO: From the ending point, you need to reach the starting point  \n',...
                ' 1. how many corners need to be added to do this (0, 1, 2, 3, 4) \n',...
                ' 2. If non-zero, enter their numbers as seen on figure-- \n']);
            noCorn=input('Enter number of corners to add:');
            sel = [];
            for zz = 1 : noCorn
                sel(zz)=input('Enter corner number):');
            end
            % if no corners, then just connect segment's end with its start.
            if isempty(sel)
                tempOuter = [segment(st:ed-1,:)
                    segment(st,:)];
            elseif sel > 0
                tempOuter = [segment(st:ed-1,:)
                    boubox(sel,:)
                    segment(st,:)];
            else
                % ignore segment
                continue;
            end
            disp('Clipping polygon with the bounding box...');
            %--Coordinates have to be positive..
            tempclippedPolygon = sutherlandHodgman(abs(tempOuter),abs(boubox));
            tempclippedPolygon(:,1) = -tempclippedPolygon(:,1);
            tempclippedPolygon = [tempclippedPolygon;tempclippedPolygon(1,:); NaN,NaN];
            
            clippedPolygon = [clippedPolygon; tempclippedPolygon];
            tempOuter = [];
            fprintf(1,' ------------------------------------------------------->\n') ;
            fprintf(1,' ------------------------------------------------------->\n') ;
        end
    end
end
