function [p,t] = General_distmesh(mapfile,bathyfile,edgelength,dist_param,...
                    wl_param,slope_param,itmax,plot_on,ini_p,fixp,poly_num)
% Function for calling grid generator for any general region in the world
% Inputs:  mapfile    : filename and path of the map created in SMS to read in
%          bathyfile  : filename of the bathymetric dataset (the version here
%                       is for a global dataset in netcdf format)
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
% Ocean Dynamics, 56(5–6), 394–415. http://doi.org/10.1007/s10236-006-0086-x
% 
    %R = 6378206.4;
    %edgelength = edgelength/R;
    %% Read the map
    [~, ~, polygon] = Read_SMS_Map( mapfile, plot_on, poly_num );
    
    %% Make new polygons based on the splitted mesh
    if ~isempty(fixp) && ~isempty(ini_p)
        % Number of individual outer polygons is equal to length I
        I = find(isnan(fixp(:,1)));
        fixp_e = fixp(1,:);
        for ii = 1:length(I)
           fixp_e = [fixp_e;fixp(I(ii)-1,:);fixp(min(I(ii)+1,end-1),:)];
        end
        fixp_e = unique(fixp_e,'rows','stable');
        idx = knnsearch(polygon.outer,fixp_e);
        ide = knnsearch(fixp,fixp_e);
        
        % Mainland goes from last fixp_e to first around anti-clockwise
        polygon.mainland = polygon.outer(idx(end):idx(1),:); 
        
        % Make outer polygon
        outer_n = [];
        for ii = 1:length(I)
           outer_n = [outer_n;fixp(ide(2*ii-1):ide(2*ii),:);...
                      polygon.outer(idx(2*ii):idx(2*ii-1),:);
                      NaN NaN];
        end
        polygon.outer = outer_n; 
        % Now remove the NaNs from fixp
        fixp(I,:) = [];
        % Remove all islands outside of outer one
        in = InPolygon(polygon.inner(:,1),polygon.inner(:,2),...
                       polygon.outer(:,1),polygon.outer(:,2));
        cn = 0;
        for ii = 1:length(in)
            if in(ii) == 1
                cn = cn + 1;
                inner_n(cn,:) = polygon.inner(ii,:);
            elseif ii > 1
                if in(ii) == 0 && in(ii-1) == 1
                   cn = cn + 1;
                   inner_n(cn,:) = [NaN, NaN];
                end
            end
        end
        if exist('inner_n','var')
            polygon.inner = inner_n;
        else
            polygon.inner = [];
        end
    end

    %% Convert all polygons to CPP
%     m_proj('Mercator' ,... %'Equidistant Cylindrical',...
%            'lon',[min(polygon.outer(:,1)) max(polygon.outer(:,1))],...
%            'lat',[min(polygon.outer(:,2)) max(polygon.outer(:,2))])
%     [ outer_x,outer_y ] = m_ll2xy( polygon.outer(:,1),polygon.outer(:,2) );
%     if ~isempty(polygon.inner)
%         [ inner_x,inner_y ] = m_ll2xy( polygon.inner(:,1),polygon.inner(:,2) );
%     else
%         inner_x = []; inner_y = [];
%     end
%     if  ~isempty(polygon.mainland)
%         [ mainland_x,mainland_y ] = m_ll2xy( polygon.mainland(:,1),polygon.mainland(:,2) );
%     else
%         mainland_x = []; mainland_y = [];
%     end
%     if ~isempty(ini_p)
%         [ini_p(:,1),ini_p(:,2)] = m_ll2xy( ini_p(:,1),ini_p(:,2) );
%     end
%     if ~isempty(fixp)
%         [fixp(:,1),fixp(:,2)] = m_ll2xy( fixp(:,1),fixp(:,2) );
%     end
    %% Make bounding box
    bounding_box = [min(polygon.outer(:,1)), min(polygon.outer(:,2)); ...
                    max(polygon.outer(:,1)), max(polygon.outer(:,2))];

    %% Make edge function interpolant
    nn = 0; % Counter for different criteria
    if slope_param > 0 || wl_param > 0 || dist_param > 0
        lon = double(ncread(bathyfile,'lon'));
        lat = double(ncread(bathyfile,'lat'));  
        bathy = double(ncread(bathyfile,'bathy'));  
        I = find(lon > min(polygon.outer(:,1)) & lon < max(polygon.outer(:,1)));
        J = find(lat > min(polygon.outer(:,2)) & lat < max(polygon.outer(:,2)));
        lon = lon(I); lat = lat(J); 
        bathy = -bathy(I,J); bathy(bathy < 0) = 0;
        [lon_g,lat_g] = ndgrid(lon,lat);
        if dist_param > 0
            nn = nn + 1;
            x_v = reshape(lon_g,[],1);
            y_v = reshape(lat_g,[],1);
            d = fd( [x_v,y_v], 1 ) ;
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
        
        % Get min of slope and wavelength and smooth
        if dist_param > 0
            hh_m = min(hh(:,:,2:end),[],3); 
        else
            hh_m = min(hh(:,:,1:end),[],3); 
        end
        
%         figure;
%         pcolor(lon,lat,hh_m')
%         shading interp

%         figure;
%         pcolor(lon,lat,hh_m1')
%         shading interp
        % Smooth by checking the adjacent length change fraction
        %[hh_m,~] = smooth2_lengthchange(hh_m,2*dist_param,1);
        %hh_m = smooth2a(hh_m,5,5);

        % Now convert hh_m into degrees from meters 
        % (estimate from 45 deg azimuth)
        [lon2,lat2,~] = m_fdist(lon_g,lat_g,45,hh_m);
        hh_m = sqrt((lon2 - lon_g).^2 + (lat2 - lat_g).^2);
        
        % Smoothing
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
        % Get min of all the criteria
        if dist_param > 0
            hh_m = min(hh_m,hh(:,:,1)); 
        end
        hh_m(hh_m < edgelength) = edgelength;
        % Make the overall interpolant
        F = griddedInterpolant(lon_g,lat_g,hh_m,'linear');
        %F = griddedInterpolant(x_g,y_g,hh_m,'linear');
        clear lon2 lat2 hh
    end

    
    %% Call distmesh        
    [p,t] = distmesh2d(@fd,@fh,@fci,...
                       edgelength,bounding_box,ini_p,fixp,itmax,plot_on);
    
    %%% Convert back to lon, lat
    %[p(:,1),p(:,2)] = m_xy2ll(p(:,1),p(:,2));
    
    if plot_on == 1
        % Plot the map
        hold on
        m_plot(polygon.outer(:,1),polygon.outer(:,2))
        m_plot(polygon.inner(:,1),polygon.inner(:,2))
    end
 
    function d = fd( p, instance ) 
        % The distance function, (needs polygons of coastline and islands)
        if instance == 1
            % When called from h func
            % First polygon is full one, second one is sans ocean boundaries
            d = dpoly( p,[polygon.outer;
                             NaN NaN;
                          polygon.inner],[polygon.mainland;
                                          polygon.inner] );
        else
            % When called from distmesh2d
            % Both polygons are full ones
            d = dpoly( p,[polygon.outer;
                             NaN NaN;
                          polygon.inner],[polygon.outer;
                                              NaN NaN;
                                          polygon.inner]); 
        end
        return ;
    end

    function h = fh( p )
        % The edge length function, desides on the edgelength at this
        % location, based on distance, wavelength and slope functions
        
        % Just call the interpolant created at beginning
        h = F(p);
        
        % Need to convert to lon, lat for interpolant
        %[p(:,1),p(:,2)] = m_xy2ll(p(:,1),p(:,2));
%         nn = 0; % Counter for different criteria
%         % Distance function
%         if dist_param > 0
%             nn = nn + 1;
%             d = fd( p, 1 ) ;
%             % Limit to dist_param*100 change between elements
%             hh(:,nn) = edgelength - dist_param*d ; 
%         end
% 
%         % Wavelength function        
%         if wl_param > 0 
%             nn = nn + 1;
%             g = 9.807;
%             period = 12.42*3600; % M2 period in seconds
%             B = Fb(p);           % Do the interpolation for depth
%             B = max(B,0);        % Make sure we have no negative depths
%             hh(:,nn) = max(edgelength,period*sqrt(g*B)/wl_param/R);
%         end
%         
%         % Slope function
%         if slope_param > 0
%             nn = nn + 1;
%             slope = Fs(p);       % Do the interpolation for slope    
%             hh(:,nn) = max(edgelength,2*pi*B./slope/slope_param/R);
%         end
           
        % Edgelength h is minimum of the considered criteria
        %h = min(hh,[],2);
        return ; 
    end

    function fci( p )
        % This function takes new points p and halves the edgefunction here
        %h = F(p);
        % Reshapre lon_g and lat_g for knnsearch
        lon_v = reshape(lon_g,[],1);
        lat_v = reshape(lat_g,[],1);
        
        % Find (2*nn_s)^2 grid points close to p
        nne = 2*nn_s;
        IDX = knnsearch([lon_v lat_v],p,'k',(2*nne)^2);
        IDX = unique(IDX(:));
        
        % Smooth the edgefunction
        hh_m1 = smooth2a(hh_m,nne,nne);
%        IDX_n = zeros(length(IDX),1); mean_v = zeros(length(IDX),1);
%         % Only keep the largest value in the nearest four, and get mean of
%         % other three
%         for i = 1:length(IDX)
%            [~,iloc] = max(hh_m(IDX(i,:))); 
%            IDX_n(i) = IDX(i,iloc); 
%            inloc = 1:4; inloc(iloc) = [];
%            mean_v(i) = mean(hh_m(IDX(i,inloc)));
%         end
%        [IDX,IA] = unique(IDX_n); mean_v = mean_v(IA);

        % Only overwrite the surrounding points
        hh_m(IDX) = max(edgelength,hh_m1(IDX));
        
        % Redo the interpolant
        F = griddedInterpolant(lon_g,lat_g,hh_m,'linear');
        
        %h1 = F(p);
        return;
    end
end