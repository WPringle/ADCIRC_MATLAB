function fid = FnGlobal_SAL_to_fort24( f24out, f14in, f15tipname, ...
                                       lonlat0, avisoloc, saldata )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolates the global SAL term onto the mesh and outputs a fort.24    %                                                                       %
%                                                                         %  
% Requires: readfort14.m                                                  %
%                                                                         % 
% Data required:                                                          %
% FES2004 loads. Source at: ftp://ftp.legos.obs-mip.fr/pub/soa/maree/...  %
%                           tide_model/global_solution/fes2004/           %
%                                                                         %
% FES2014 loads. Source at: ftp://ftp.legos.obs-mip.fr/pub/...            %
%                           FES2012-project/data/LSA/FES2014/             %
%                                                                         %
% Created by William Pringle Oct 20 2016 for FES2004 SAL                  %
% Updated by William Pringle Oct 28 2016 for FES2014 SAL                  %
% Updated by Dam Wa for make into a function                              %
%                                                                         %
% Run example:                                                            %
% f14 = 'fort.14' ;                                                       %
% f24 = 'fort.24'                                                         %
% latlon0 = [75.214667 -31.172085] ;                                      %
% saldat = 'FES2014' ;                                                    %
% avisoloc = './AVISO_DIREC'                                              %
% f15tipname = { 'M2', 'O1', 'S2', 'N2', 'K2', 'K1', 'Q1', 'P1'}          %
%                                                                         %
% FnGlobal_SAL_to_fort24(f24, f14, f15tipname, lonlat0, avisoloc, saldata)%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ll0 = lonlat0(1) ;
if ( ll0 < 0 ) 
    ll0 = ll0 + 360 ; 
end
lon0 = ll0*pi/180 ; lat0 = lonlat0(2)*pi/180 ; 

R = 6378206.4; % earth radius
              
% Constituents that we want in the order that we want
const = {'M2','S2','K1','O1','N2','K2','P1','Q1','M4'};
frequency = [0.000140518902509,0.000145444104333,0.000072921158358,...
             0.000067597744151,0.000137879699487,0.000145842317201,...
             0.000072522945975,0.000064958541129,0.000281037805018];
alpha = {'M2 SAL','S2 SAL','K1 SAL','O1 SAL',...
         'N2 SAL','K2 SAL','P1 SAL','Q1 SAL','M4 SAL'};


ntip = length(f15tipname) ;
iconsort = zeros(ntip,1) ;
for icon = 1: ntip
    strcon = strtrim(char(f15tipname{icon})) ;
    
    for ic = 1: length(const)
        if ( strcmpi(strcon,char(const{ic})) )
            iconsort(icon) = ic ;
            break ; 
        end
    end
end
iconsort    

% choose tidal database file names and directories
%database = 'FES2004';
% database = 'FES2014';
database = strtrim(upper(saldata)) ;
% direc    = 'E:\Global_Data\AVISO_TIDES\';
direc = strtrim(avisoloc) ;

% input fort.14 name
% fort14    = '../IDIOMS_v5.18.grd';
fort14 = strtrim(f14in) ; 

% output fort.24 name
% fort24    = ['fort.24.' database];
fort24 = [strtrim(f24out) '.' database] ;

% % Load tide grid data 
if strcmp(database,'FES2004')
    tide_grid     = [direc '/' database '/SAL/load.k1.nc'];
    tide_prefix   = [direc '/' database '/SAL/load.'];
    tide_suffix   = '.nc';
    if ( ispc )
        tide_grid     = [direc '\' database '\SAL\load.k1.nc'];
        tide_prefix   = [direc '\' database '\SAL\load.'];
        % tide_suffix   = '.nc';
    end
    
    %ncdisp(tide_grid);
    lon = ncread(tide_grid,'lon');
    lat = ncread(tide_grid,'lat');
elseif  strcmp(database,'FES2014')
    tide_grid     = [direc '/' database '/SAL/K1_sal.nc'];
    tide_prefix   = [direc '/' database '/SAL/'];
    tide_suffix   = '_sal.nc';
    if ( ispc )
        tide_grid     = [direc '\' database '\SAL\K1_sal.nc'];
        tide_prefix   = [direc '\' database '\SAL\'];
        tide_suffix   = '_sal.nc';
    end
    
    %ncdisp(tide_grid);
    lon = ncread(tide_grid,'longitude');
    lat = ncread(tide_grid,'latitude');
    [lon,lat] = ndgrid(lon,flipud(lat));
end
% CPP Conversion of lat/lon
lon = lon * pi/180; lat = lat * pi/180; 
x = R * (lon - lon0) * cos(lat0);
y = R * lat;

% % Get mesh info
if strcmp(fort14(end-2:end),'grd') || strcmp(fort14(end-1:end),'14')
    [~,VX,~,~,~,~] = readfort14( fort14 , 0) ;
else strcmp(fort14(end-2:end),'mat')
    load(fort14)
end 
% Ensuring 0 to 360
VX(VX(:,1) < 0,1) = VX(VX(:,1) < 0,1) + 360;
% Doing the CPP conversion
xx = VX(:,1) * pi/180; yy = VX(:,2) * pi/180; 
xx = R * (xx - lon0) * cos(lat0);
yy = R * yy;

% % Now interpolate onto grid and write out to fort.24 type file

fid = fopen(fort24,'w');

nnodes = length(VX) ;
kvec = [1:nnodes]' ; 
for icon = 1: ntip
    j = iconsort(icon) ; 
% for j = 1:length(const)

    % The current consituent filename
    if strcmp(database,'FES2004')
        tide = [tide_prefix lower(const{j}) tide_suffix];
        % Get amp and phase
        Ha = ncread(tide,'Ha');
        Hg = ncread(tide,'Hg');
    elseif  strcmp(database,'FES2014')
        tide = [tide_prefix const{j} tide_suffix];
        % Get amp and phase
        Ha = ncread(tide,'SAL_amplitude');
        Hg = ncread(tide,'SAL_phase');
        Ha = fliplr(Ha);
        Hg = fliplr(Hg);
    end

%     figure;
%     [c,h] = contour(lon,lat,Ha);
%     clabel(c,h)
%     colorbar
%     figure;
%     [c,h] = contour(lon,lat,Hg,0:30:360);
%     clabel(c,h)
%     colorbar
    
    Hg(Hg > 180) = Hg(Hg > 180) - 360; % move to -180 - 180
    Hg = Hg*pi/180; %radians
    
    % Convert to complex number for interpolation
    z = Ha.*exp(Hg*1i);
    
    % Do the gridded Interpolation
    F = griddedInterpolant(x,y,z,'linear');
    Z = F(xx,yy);  

    % Convert back to amp and phase
    amp = abs(Z);
    phs = angle(Z)*180/pi;
    
    % Convert to 0 to 360;
    % phs_b that is positive 0 - 180 stays same.
    % phs_b that is negative -180 - 0 becomes 180 - 360
    phs(phs < 0) = phs(phs < 0) + 360;
    
    % Print out interpolated results
    
    fprintf('Wrting SAL %s data \n', char(const{j})) ; 
    % The constituent details
    fprintf(fid,'%s \n',cell2mat(alpha(j))) ;
    fprintf(fid,'%17.15f \n',frequency(j)) ;
    fprintf(fid,'%d \n',1) ;  
    fprintf(fid,'%s \n',cell2mat(const(j))) ;
    
    % Loop over the nodes of the mesh
    % for k = 1:length(amp)
    %    fprintf(fid,'%d \t %12.6f  %12.6f \n',k,amp(k),phs(k));
    % end
    
    fprintf(fid,'%d \t %12.6f  %12.6f \n',[kvec amp phs]');

end
fclose(fid);
end
