function Nb_av = Compute_Nb_greens_av(VX,B,Nm,Nb,omega,proj,hcut,acut)
% Computes Gaussian smoothing of Nb scaled by cutoff length scale used in 
% Nycander (2005) formula
% 
% Inputs:  VX   - the vertex positions (N x 2)
%          B    - the bathymetric depth at each node (N x 1)
%          Nm   - the depth-averaged Brunt–Väisälä buoyancy frequency at
%                 each node (N x 1)
%          Nb   - the Brunt–Väisälä buoyancy frequency at the seabed at
%                 each node (N x 1)
%         omega - the circular frequency (rad/s) of the tidal 
%                 component of interest (usually the M2 one)
%          proj - a string defining the projection to use, e.g. 'Mercator'.
%                 See m_proj('get') for other options in m_map 
%         hcut  - The cutoff height (only smooth Nb.. when depth > hcut)
%         acut  - The cutoff radius (only smooth Nb.. when a > acut)
%
% Outputs: Nb_av  - Gaussian smoothed Nb (N x 1)
%
% Authors : William Pringle Oct 4 2017
%
% Requires : - m_map toolbox for projection

% Reference: Nycander, J. (2005). Generation of internal waves in the deep
%            ocean by tides. Journal of Geophysical Research, 
%            110(C10), C10028. http://doi.org/10.1029/2004JC002487

%% Set some constants
% Earth Radius
R = 6378206.4;
% Coriolis coefficient
psi = 2*7.29212d-5;
% beta coefficient
beta = 1.455;

% Some arbitrary settings for the numerics 
rcutfac = 2.5 ; % Elements within |r - r'| < rcutfac*a will 
                % be used in the smoothing of Nb
    
%% Computing the cutoff length of green's function
f = psi*sind(VX(:,2)) ; % Coriolis coefficients
dom = pi*sqrt(omega^2 - f.^2);
a = (beta./dom).*B.*Nm ;      % the cutoff length

%% Doing the projection
% Setting up projection
m_proj(proj,'lon',[ min(VX(:,1)) max(VX(:,1))],...
            'lat',[ min(VX(:,2)) max(VX(:,2))])
% Doing conversion
[XX(:,1),XX(:,2)] = m_ll2xy(VX(:,1),VX(:,2));      
XX = R * XX; % Convert to actual distances 
np = length(XX) ;  % Number of points

% Initialising the arrays
if length(hcut) > 1
    idv = find( (B > hcut(1)) & (B < hcut(2)) & (a > acut) ) ;
else
    idv = find( (B > hcut)  & (a > acut) ) ;
end
    
nj = length(idv) ;

Nb_av  = zeros(np,1)  ;

%% Doing the computation
NbTemp   = zeros(nj,1) ; 

% Gaussian function
gaussC = @(dist, sigma2) pi^(-3/2) * (exp(-dist.^2./sigma2)); 
%greensF = @(dist) 1./dist - (sqrt(pi)/2) * ...
%                  (exp(-dist.^2./8)) .* besseli(0,dist.^2./8);


% set up kdtree
Mdl = KDTreeSearcher(XX); 
% sort by radius
[aiv, ia] = sort(a(idv));
xiv = XX(idv(ia),:);
% number of times to do rangesearch
blocks = ceil(nj*rcutfac^2*1.5e-6);
nn = 0;
disp(['number of blocks = ' num2str(blocks)])
for ii = 1:blocks
    tic
    ns = int64((ii-1)*nj/blocks)+1;
    ne = int64(ii*nj/blocks);
    % Do range search
    [Idx, Dist] = rangesearch(Mdl,xiv(ns:ne,:),rcutfac*aiv(ne));  
    % Loop through Idx
    for jj = 1:length(Idx)
        nn = nn + 1;
        val = gaussC(Dist{jj}(Dist{jj} < rcutfac*aiv(nn)), aiv(nn)^2); % Get the weights
        %val = val/norm(val,1); % divide by the norm
        %val2 = greensF(Dist{jj}(Dist{jj} < rcutfac*aiv(nn))/aiv(nn)); % Get the weights
        NbTemp(nn) = val*Nb(Idx{jj}(Dist{jj} < rcutfac*aiv(nn))); %sum the Nb with weights
        % NbTemp(nn) = mean(Nb(Idx{jj}(Dist{jj} < aiv(nn))));
    end
    disp(['finished block no. = ' num2str(ii)])
    toc
end
% Enter into global arrays
Nb_av(idv(ia)) = NbTemp ;  
%EOF
end