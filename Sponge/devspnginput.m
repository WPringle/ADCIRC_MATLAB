% This script calls functions for making the sponge zone, outputting the
% lat, lon in that sponge zone, outputting the fort13 of sigma coefficients
% and finally making the fort.53001 and fort.54001 of TPXO eta, u, & v
% 
% Note that we use uH from TPXO and then calculate u using our bathy
%
clearvars; clc;

% Extract tpxo folder
hmtpxo = './OTPSEXTRACTDAT/' ;

% set suffix of file names
hsuf = ['.dat'] ; 
usuf = ['U.dat'] ;
vsuf = ['V.dat'] ; 

% our grid
f14file = '../INDPAC_SP.grd' ;

% sponge boundaries (select which boundaries to set a sponge at
sbnd = [1,3];

% set sponge type 
spngtype = 'poly'; 
F = 20; % parameters 
rat = 1/3; % sponge type
period = 12.42*3600; %period of M2 wave
frac   = 0.1; % length of the sponge zone/length of M2 wave

% mindepth allowable for dividing TPXO uH by our bathymetry H
dmin = 10;

% NB: First iteration set write == 1, get the spng_lat_lon and fort.13
% file. Interpolate OTPS onto spng_lat_lon elsewhere then do second
% iteration setting write = 0.
write = 0; % write the spong_lat_lon or not

%% Make sponge and output fort 13
f13file = ['fort.13.sp_F' num2str(F) '_R' num2str(rat)];

% make the fort13 with the sponge coefficients
[sponge,opedat,boudat,pv,B] = makefort13sponge(f13file,f14file,period,...
                                           frac,spngtype,F,rat,sbnd,write);
if write == 1
    return;
end

%% Make the fort.15, and get info from tpxO8
% Get the constituents info from a fort.15
f15dat = readfort15_to_NBFR( 'fort.15' ) ;

% loop over all the constituents, output the TPXO8 boundary data and the
% fort.53001 and fort.54001 sponge data

% Loop over constituents and get data
Ha  = cell(f15dat.nbfr,1);
U   = cell(f15dat.nbfr,1);
V   = cell(f15dat.nbfr,1);
for k = 1:f15dat.nbfr
    % Get the amp
%     T = readtable([hmtpxo strtrim(f15dat.bountag(k).name) hsuf],...
%                   'HeaderLines',2);
    [T,~] = readotpsout( [hmtpxo ...
                           lower(strtrim(f15dat.bountag(k).name)) hsuf] ); 
              
    % Get nearest interp for -99999 values and save for fort.53001
    Ha{k} = nearestinterpOTPS(T);
%     J = find(table2array(T(:,3)) == -99999 ); 
%     if ~isempty(J)
%         K = find(table2array(T(:,3)) ~= -99999 );
%         IDX = knnsearch([T.Lon(K),T.Lat(K)],[T.Lon(J),T.Lat(J)]);
%         T(J,3) = T(K(IDX),3);
%         T(J,4) = T(K(IDX),4);     
%     end
    % Save for fort.53001
    %Ha{k} = table2array([T(:,3), T(:,4)]);
    
    % Enter into the f15 structure
    f15dat.opeemoefa(k).val = [];
    for op = 1:opedat.nope
        I = knnsearch([T(2,:)' T(1,:)'],pv(opedat.nbdv(1:opedat.nvdll(op),op),:));
        f15dat.opeemoefa(k).val = [f15dat.opeemoefa(k).val; T(3:4,:)'];
    end
    
    % Get the U
    [T,~] = readotpsout( [hmtpxo ...
                           lower(strtrim(f15dat.bountag(k).name)) usuf] );
    % Get nearest interp for -99999 values and save for fort.54001
    U{k} = nearestinterpOTPS(T);
    % Get the V
    [T,~] = readotpsout( [hmtpxo ...
                           lower(strtrim(f15dat.bountag(k).name)) vsuf] );
    % Get nearest interp for -99999 values and save for fort.54001
    V{k} = nearestinterpOTPS(T);
end
% Output the open boundary stuff to fort15
writefort15_to_BND( 'fort.15.tpxO',f15dat ) ;

%% Make fort.53001 and fort.54001
% open files
f53 =  fopen('fort.53001','w') ;
f54 =  fopen('fort.54001','w') ;
% Enter header 
% NFREQ
fprintf(f53,'%d \n', f15dat.nbfr ) ;
fprintf(f54,'%d \n', f15dat.nbfr) ;
for k = 1: f15dat.nbfr
    fprintf(f53,'%14.9f %f %f %s\n', f15dat.bounspec(k).val, f15dat.bountag(k).name ) ; 
    fprintf(f54,'%14.9f %f %f %s\n', f15dat.bounspec(k).val, f15dat.bountag(k).name ) ; 
end
% NP
ne = length(sponge(1).idx) + length(sponge(3).idx);
fprintf(f53,'%d \n', ne) ; 
fprintf(f54,'%d \n', ne) ; 
nn = 0;
for op = 1:opedat.nope
    for i = 1:length(sponge(op).idx)
        nn = nn + 1;
        % Skipping open boundaries where we don't have sponge
        if ~any(op == sbnd); continue; end
        % print node number
        fprintf(f53, '%d \n', sponge(op).idx(i)) ; 
        fprintf(f54, '%d \n', sponge(op).idx(i)) ; 

        H = B(sponge(op).idx(i));
        % print amp and phases
        for k = 1: f15dat.nbfr
            if H <= 0
                % Set to zero
                fprintf(f53, '%14.9e %14.9e \n', 0,0) ;
                fprintf(f54, '%14.9e %14.9e %14.9e %14.9e \n', ...
                        0,0,0,0);
            else
                fprintf(f53, '%14.9e %14.9e \n', Ha{k}(nn,:)) ;
                % Divide uH from TPXO8 by H for the amp only
                fprintf(f54, '%14.9e %14.9e %14.9e %14.9e \n', ...
                        U{k}(nn,:)./[max(dmin,10) 1],...
                        V{k}(nn,:)./[max(dmin,10) 1]);
            end
        end
    end
end
fclose(f53);
fclose(f54);


           

