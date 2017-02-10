function [p,t] = General_distmesh(mapfile,bathyfile,edgelength,dist_param,...
                                  wl_param,slope_param,itmax,plot_on)
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
    R = 6378206.4;
    edgelength = edgelength/R;
    % Read the map
    [~, ~, polygon] = Read_SMS_Map( mapfile, plot_on );

    % Convert all polygons to CPP
    m_proj('Equidistant Cylindrical',...
           'lon',[min(polygon.outer(:,1)) max(polygon.outer(:,1))],...
           'lat',[min(polygon.outer(:,2)) max(polygon.outer(:,2))])
    [ outer_x,outer_y ] = m_ll2xy( polygon.outer(:,1),polygon.outer(:,2) );
    [ inner_x,inner_y ] = m_ll2xy( polygon.inner(:,1),polygon.inner(:,2) );
    [ mainland_x,mainland_y ] = m_ll2xy( polygon.mainland(:,1),polygon.mainland(:,2) );
    
    % Make bounding box
    bounding_box = [min(outer_x), min(outer_y); ...
                    max(outer_x), max(outer_y)];

    % Load bathy and make interpolants
    if slope_param > 0 || wl_param > 0
        lon = double(ncread(bathyfile,'lon'));
        lat = double(ncread(bathyfile,'lat'));  
        bathy = double(ncread(bathyfile,'bathy'));  
        I = find(lon > min(polygon.outer(:,1)) & lon < max(polygon.outer(:,1)));
        J = find(lat > min(polygon.outer(:,2)) & lat < max(polygon.outer(:,2)));
        lon = lon(I); lat = lat(J); bathy = bathy(I,J);
        [lon_g,lat_g] = ndgrid(lon,lat);
        [ x_g,y_g ] = m_ll2xy( lon_g,lat_g );
    
        % Making the bathymetric interpolant
        Fb = griddedInterpolant(x_g,y_g,-bathy,'linear');
    end
        
    % Making the slope interpolant
    if slope_param > 0 
        dx = R*(x_g(2,1)-x_g(1,1)); dy = R*(y_g(1,2) - y_g(1,1));
        [b_x,b_y] = gradient(bathy,dx,dy);
        b_slope = sqrt(b_x.^2 + b_y.^2);
        Fs = griddedInterpolant(x_g,y_g,b_slope,'nearest');
    end
    
    % Call distmesh        
    [p,t] = distmesh2d(@fd,@fh,edgelength,bounding_box,[],itmax,plot_on) ;
    
    % Convert back to lon, lat
    [lon_p,lat_p] = m_xy2ll(p(:,1),p(:,2));
    p = [lon_p,lat_p];
    
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
            d = dpoly( p,[outer_x,outer_y;
                             NaN NaN;
                          inner_x,inner_y],[mainland_x,mainland_y;
                                            inner_x,inner_y] );
        else
            % When called from distmesh2d
            % Both polygons are full ones
            d = dpoly( p,[outer_x,outer_y;
                             NaN NaN;
                          inner_x,inner_y],[outer_x,outer_y;
                                              NaN NaN;
                                            inner_x,inner_y]); 
        end
        return ;
    end

    function h = fh( p )
        % The edge length function, desides on the edgelength at this
        % location, based on distance, wavelength and slope functions
        
        nn = 0; % Counter for different criteria
        % Distance function
        if dist_param > 0
            nn = nn + 1;
            d = fd( p, 1 ) ;
            % Limit to dist_param*100 change between elements
            hh(:,nn) = edgelength - dist_param*d ; 
        end

        % Wavelength function        
        if wl_param > 0 
            nn = nn + 1;
            g = 9.807;
            period = 12.42*3600; % M2 period in seconds
            B = Fb(p);           % Do the interpolation for depth
            B = max(B,0);        % Make sure we have no negative depths
            hh(:,nn) = max(edgelength,period*sqrt(g*B)/wl_param/R);
        end
        
        % Slope function
        if slope_param > 0
            nn = nn + 1;
            slope = Fs(p);       % Do the interpolation for slope    
            hh(:,nn) = max(edgelength,2*pi*B./slope/slope_param/R);
        end
           
        % Edgelength h is minimum of the considered criteria
        h = min(hh,[],2);
        return ; 
    end
end