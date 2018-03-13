function obj = Calc_IT_Fric(obj,varargin)
% f13dat = Calc_IT_Fric(obj,N_filename,type,C_it,MinDepth,crit)
% Input a msh class object with bathy and slope data, get the values of N
% over the depth based on N_filename (if is not empty) and calculate the
% internal tide friction based on the input parameters   
% 
%  Inputs:      1) .mat files of N values at constant contours    %
%               2) Unstructured grid mesh with bathymetry         %                           
%  Outputs:     A fort.13 formatted file for use in ADCIRC/SMS    %
%  Project:     Indian Ocean and Marginal Seas                    %
%  Author:      William Pringle                                   %
%  Created:     Oct 5 2016                                        %
%  Updated:     Oct 24 2016, Dec 14 2016, Feb 25 2017             %
%  Requires:    functions - readfort14, Compute_Nb_Nm, m_proj,    %
%               m_ll2xy, ADCIRC_Bath_Slope, Compute_J_Nycander    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Some constants
% Choose projection type (can be any, not restricted to evenly-gridded data)
proj = 'Mercator';
% Radius of earth for conversion to actual distances
R = 6378206.4; %[m]
% Coriolis coefficient
psi = 2*7.29212d-5;
% Choose tidal frequency to compute for (in rads^{-1}) - usually M2 freq.
omega = 2*pi/(12.4206012*3600); % M2 tidal frequency

%% Test optional arguments
% default
type = 'directional';
crit = 'Nb';
C_it = 0.25;
MinDepth = 100; % m 
N_filename = []; % no N data (just print out slopes multiplied by C_it)
if ~isempty(varargin)
    names = {'type','cutoff_depth','Cit','Nfname'};
    for ii = 1:length(names)
        ind = find(~cellfun(@isempty,strfind(varargin(1:2:end),names{ii})));
        if ~isempty(ind)
            if ii == 1
                type = varargin{ind*2}; 
            elseif ii == 2
                MinDepth = varargin{ind*2};
            elseif ii == 3
                C_it = varargin{ind*2};
            elseif ii == 4 
                N_filename = varargin{ind*2};
            end
        end    
    end
end

% Coriolis
f = psi*sind(obj.p(:,2));           
    
%% Load the constant contours of N values and compute Nb and Nmean
if ~isempty(N_filename)
    load(N_filename);  
    [Nb,Nm,Nmw] = Compute_Nb_Nm_Gridded(obj.p(:,1),obj(:,2),obj.b,z,N,lon,lat);                   
end

%% Getting the J stuff if required (tensor type)
if strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
        % Compute gradients of J from Nb, Nm and bathymetry B, and grid  
        [J,dJ] = Compute_J_Nycander(obj.t,obj.p,obj.b,Nm,omega,...
                                       2,MinDepth,proj,[],4);                    
end
%
H2_mesh = obj.bx.^2 + obj.by.^2;

%% Calculate F_it from Nb, Nm, Jx, Jy, and H2_mesh or as reqd.
if isempty(N_filename)
    
else
    Nb_t = Nb;
    if strcmp(crit,'Nmw')
        Nb_t = Nmw;
    end
    if strcmp(type,'scalar')
       F_it = C_it * sqrt((Nb_t.^2 - omega^2).*(Nm.^2 - omega^2)).*H2_mesh/omega;
    elseif strcmp(type,'directional')
       F_it = C_it * sqrt((Nb_t.^2 - omega^2).*(Nm.^2 - omega^2))/omega;
    elseif strcmp(type,'tensor') || strcmp(type,'tensor_to_scalar')
       F_it = C_it * Nb_t/(4*pi).*sqrt(1-f.^2/omega^2)./obj.b;
    end
    F_it = real(F_it); % in case becomes complex due to minus square root

    %% Compute criticality and normalise the friction
    alpha2  = (omega^2 - f.^2)./(Nb_t.^2 - omega^2);
    if ~strcmp(crit,'none')
        gamma2 = max(H2_mesh./alpha2,1);
        % Normalise F_it by criticality
        F_it = F_it./gamma2;
        % Make sure that if alpha2 < 0 that F_it is 
        % set equal to 0 since real part would be 0.
        F_it(alpha2 < 0) = 0;
    end
end

% Cut it off at the MinDepth
F_it(obj.b < MinDepth) = 0;

% Make zero for no slopes
F_it(obj.bx == 0 & obj.bx == 0) = 0;
if strcmp(type,'tensor')
    F_it(dJ(:,1) == 0 & dJ(:,1) == 0) = 0;
end
%
%% Make into f13 struct
if isempty(obj.f13)
    
end
% Print header
fprintf(fid,'%s \n','Spatial attributes description') ;
fprintf(fid,'%d \n',length(F_it)) ;   
fprintf(fid,'%d \n',1) ;  
fprintf(fid,'%s \n','internal_tide_friction') ;
fprintf(fid,'%s \n','1/time') ;
if strcmp(type,'tensor') || strcmp(type,'directional')
    fprintf(fid,'%d \n',3) ;
    fprintf(fid,'%f %f %f\n',0.0,0.0,0.0) ;  
else
    fprintf(fid,'%d \n',1) ;  
    fprintf(fid,'%f \n',0.0) ;  
end
fprintf(fid,'%s \n','internal_tide_friction') ;
%
% Number of nodes not default (0.0)
ne = length(find(F_it > 0));
fprintf(fid,'%d \n',ne) ; 
% Print out list of nodes for each
for k = 1:length(F_it)
    if F_it(k) > 0
        if strcmp(type,'tensor')
            % Output the C_11, C_12 = C_21, C_22 for the tensor
            C_11 = 2*F_it(k)*dJx(k)*Hy(k);
            C_22 = 2*F_it(k)*dJy(k)*Hx(k);
            C_21 = F_it(k)*(dJx(k)*Hy(k) + dJy(k)*Hx(k));
            %C_11 = 2*F_it(k)*dJ(k,1)*dh(k,1);
            %C_22 = 2*F_it(k)*dJ(k,2)*dh(k,2);
            %C_21 = F_it(k)*(dJ(k,1)*dh(k,2) + dJ(k,2)*dh(k,1));
            fprintf(fid,'%d \t %15.9e %15.9e %15.9e \n',...
                    k,C_11,C_22,C_21);
        elseif strcmp(type,'directional')
            % Output the C_11, C_12 = C_21, C_22 for the tensor
            C_11 = F_it(k)*Hx(k)^2;
            C_22 = F_it(k)*Hy(k)^2;
            C_21 = F_it(k)*Hx(k)*Hy(k);
            fprintf(fid,'%d \t %15.9e %15.9e %15.9e \n',...
                    k,C_11,C_22,C_21); 
        else
            % Only need the F_it component
            fprintf(fid,'%d \t %15.9e \n',k,F_it(k));
            %fprintf(fid,'%d \t %15.9e \n',k,Nb(k));
        end
    end
end
fclose(fid);


