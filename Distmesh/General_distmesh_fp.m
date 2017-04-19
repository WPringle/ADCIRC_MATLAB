function [p,t] = General_distmesh_fp(meshfile,contourfile,bathyfile,bbox,...
      edgelength,dist_param,wl_param,slope_param,ub,minL,itmax,plot_on)
% Function for calling grid generator for any general region in the world
% Inputs:  meshfile   : filename of the ADCIRC mesh created for the
%                       oceanside
%          contourfile: filename of the shapefile of the upper bound contour
%          bathyfile  : filename of the bathymetric dataset (the version
%                       here is for one created from a DEM in gdal)
%          edglength  : The minimum (and initial) edgelength of the mesh
%          dist_param : The percent that the edgelength should change 
%                       with distance. Set 0 to ignore
%          wl_param   : Parameter in wavelength function, (set 0 to ignore)
%                       edgelength = T/wl_param*sqrt(gH), where T is the
%                       period of the M2 wave and H is the bathymetric depth
%          slope_param: Parameter in slope function, (set 0 to ignore)
%                       edgelength = 2*pi/slope_param*H/abs(grad(H)), 
%                       where H is the bathymetric depth
%          itmax      : Maximum number of iterations in distmesh function
%          plot_on    : choose whether to plot map and mesh during calc or
%                       not (Yes = 1, No = 0)
%          ini_p      : optional initial distribution of p
%          fix_p      : optional fixed node points
%          poly_num   : the polygon number inside the mapfile
%
% Outputs: p          : the node positions in lon,lat
%          t          : the triangulation of the mesh
%
% By William Pringle, Damrongsak Wirasaet Feb 2017
%
% Reference for wavelength and slope function:
% Lyard, F., Lefevre, F., Letellier, T., & Francis, O. (2006). 
% Modelling the global ocean tides: modern insights from FES2004. 
% Ocean Dynamics, 56(5�6), 394�415. http://doi.org/10.1007/s10236-006-0086-x
% 
    %% Read the mesh
    [ev,pv,~,opendat,boudat,~] = readfort14(meshfile);
    % initialise segment and polygon
    segment = []; polygon = [];
    [etbv,~,~,~] = extdom_edges( ev, pv ) ;% the edges in the external domain
    iedbeg  = knnsearch(pv(etbv(1,:),:),pv(opendat.nbdv(1),:)); % find which edge to start
    polygon = extdom_polygon( etbv, pv, iedbeg, 1, 2 ) ; % find the outer polygon defining the domain
    polygon = [polygon; NaN NaN];
    % Make boudat into segment and polygons based on min criteria
    for i = 1:boudat.nbou
        nodes = boudat.nbvv(1:boudat.nvell(i),i);
        if boudat.nvell(i) > minL
            % Add to segment 
            segment = [segment;NaN NaN; pv(nodes,:)];
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
    polygon = [polygon;NaN NaN];
    % make fixp = to the segment
    fixp  = segment; 
    % get rid of the nans
    fixp=fixp(~any(isnan(fixp),2),:);
    ini_p = [];

    %% Read floodplain contour
    % add floodplain contour to segment   
    S = shaperead(contourfile);
    % Get rid of unwanted components
    F = fieldnames(S);
    D = struct2cell(S);
    S = cell2struct(D(3:4,:), F(3:4));
    for s = 1:length(S)
        x_n = S(s).X; y_n = S(s).Y; 
        x_n = x_n(~isnan(x_n)); y_n = y_n(~isnan(y_n));
        m_d = mean(abs(diff([x_n; y_n],[],2)),2); m_d = norm(m_d,2);
        % Ignore small length shapes
        if length(x_n) < edgelength/m_d*7; continue; end
        I = find(x_n > bbox(1,1) & x_n < bbox(1,2) & ...
                 y_n > bbox(2,1) & y_n < bbox(2,2));
        % Ignore small length shapes within bbox
        if length(I) < edgelength/m_d*7; continue; end        
        segment = [segment; NaN NaN; x_n(I)' y_n(I)'];
    end
    % Build kdtree with segment points for distance calc algorithm 
    mdl = KDTreeSearcher(segment); 

    %% Read DEM and get bathy interpolant for fd
    lon = double(ncread(bathyfile,'lon'));
    lat = double(ncread(bathyfile,'lat'));  
    bathy = double(ncread(bathyfile,'Band1'));  %Band1 if created in gdal
    I = find(lon > bbox(1,1) & lon < bbox(1,2));
    J = find(lat > bbox(2,1) & lat < bbox(2,2));
    lon = lon(I); lat = lat(J); bathy = bathy(I,J);
    [lon_g,lat_g] = ndgrid(lon,lat);
    Fb = griddedInterpolant(lon_g,lat_g,bathy);
    
    %% Make edge function interpolant
    nn = 0; % Counter for different criteria
    if dist_param > 0
        nn = nn + 1;
        x_v = reshape(lon_g,[],1);
        y_v = reshape(lat_g,[],1);
        d = fd( [x_v,y_v], 0 ) ;
        % reshape back
        d = reshape(d,length(I),[]);
        % add into edge function
        hh(:,:,nn) = edgelength - dist_param*d ; 
        clear x_v y_v d
    end    
    if wl_param > 0
        nn = nn + 1;
        g = 9.807;
        period = 12.42*3600; % M2 period in seconds
        % wavelength func
        hh(:,:,nn) = period*sqrt(g*bathy)/wl_param;
    end
    % Making the bathymetric interpolant
    %Fb = griddedInterpolant(x_g,y_g,-bathy,'linear');
    if slope_param > 0 
        nn = nn + 1;
        lon2 = [lon_g(2:end,:); lon_g(end,:)+(lon(end)-lon(end-1))];
        lat2 = [lat_g(:,2:end) lat_g(:,end)+(lat(end)-lat(end-1))];
        [dx,~,~] = m_idist(lon_g,lat_g,lon2,lat_g);
        [dy,~,~] = m_idist(lon_g,lat_g,lon_g,lat2);
        % Now get the coordinate matrices for dx and dy
        dx = cumsum(dx,1); dx = cumsum(dx,2);
        dy = cumsum(dy,1); dy = cumsum(dy,2);
        [b_y,b_x] = gradient(bathy,dy,dx);
        b_slope = sqrt(b_x.^2 + b_y.^2);
%             figure;
%             pcolor(lon,lat,b_slope')
%             shading interp
%             caxis([0 0.1])
        % Smooth slope using three point Shapiro filter fnum times
        fnum = 1;
        % add buffer
        for n = 1:fnum
            b_slope = [b_slope(:,1) b_slope b_slope(:,end)];
            b_slope = [b_slope(1,:); b_slope; b_slope(end,:)];
        end
        % perform filter
        for n = 1:fnum
            b_slope = 0.25*(b_slope(:,1:end-2) + 2*b_slope(:,2:end-1) + b_slope(:,3:end));
            b_slope = 0.25*(b_slope(1:end-2,:) + 2*b_slope(2:end-1,:) + b_slope(3:end,:));
        end
        % slope func
        hh(:,:,nn) = 2*pi*bathy./b_slope/slope_param;
        clear lon2 lat2 dx dy b_y b_x b_slope
    end

    % Get min of slope and wavelength, setting equal to hh_m
    if wl_param > 0 || slope_param > 0
        if dist_param > 0 
            hh_m = min(hh(:,:,2:end),[],3); 
        else
            hh_m = min(hh(:,:,1:end),[],3); 
        end
        % Now convert hh_m into degrees from meters 
        % (estimate from 45 deg azimuth)
        [lon2,lat2,~] = m_fdist(lon_g,lat_g,45,hh_m);
        lon2(lon2 > 180) = lon2(lon2 > 180) - 360;
        hh_m = sqrt((lon2 - lon_g).^2 + (lat2 - lat_g).^2);
        clear lon2 lat2
        % Smoothing if slope_param exists
        if slope_param > 0 
            %
            %OPTIONS.Weight = 'cauchy';
            %OPTIONS.MaxIter = 100;
            %[hh_m1,S,EXITFLAG] = smoothn(hh_m,'robust',OPTIONS);
            %nx = 2*round((size(hh_m,2)*0.005+1)/2)-1;
            %nz = 2*round((size(hh_m,1)*0.005+1)/2)-1;
            %hh_m2 = savitzkyGolay2D_rle_coupling(size(hh_m,2),size(hh_m,1),...
            %                                     hh_m,nx,nz,2);
            I = bad_length_change(hh_m,dist_param*2); nn_s = 1;
            while length(I)/length(hh_m(:)) > 0.005
                hh_m1 = smooth2a(hh_m,nn_s,nn_s);
                nn_s = nn_s + 1;
                I = bad_length_change(hh_m1,dist_param*2);
                if nn_s > max(size(hh_m))*0.005
                    break;
                end
            end
            if exist('hh_m1','var')
                hh_m = hh_m1;
            end
        end
        % Get min of all the criteria
        if dist_param > 0
            hh_m = min(hh_m,hh(:,:,1)); 
        end
    else
        hh_m = squeeze(hh);
    end
    hh_m(hh_m < edgelength) = edgelength;    
    % Make the overall interpolant
    F = griddedInterpolant(lon_g,lat_g,hh_m,'linear');
    clear hh
    
    %% Call distmesh        
    [p,t] = distmesh2d(@fd,@fh,...   
                       edgelength,bbox',ini_p,fixp,itmax,plot_on);
    if plot_on == 1
        % Plot the map
        hold on
        m_plot(polygon.outer(:,1),polygon.outer(:,2))
        m_plot(polygon.inner(:,1),polygon.inner(:,2))
    end
 
    function d = fd( p, instance )
        % The distance function 
        % needs segment to get distance
        % smallpolygons to make sure nodes inside them are discarded
        % Fb and bounds to test whether depth in desired range
        % NOTE: Pass mdl here instead of segment
        d = dpoly_fp2(p,mdl,polygon,Fb,ub);
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
