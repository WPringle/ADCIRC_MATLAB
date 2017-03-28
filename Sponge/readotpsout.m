function [otps,vardes] = readotpsout( otpsfile )
% read OTPS outputfile

%otpsfile = 'k1.out' ;

fid = fopen( otpsfile ) ; 

% ignore first three line
line = fgetl( fid ) ;
line = fgetl( fid ) ;
vardes = fgetl( fid ) ;

otps = fscanf(fid, '%f %f %f %f\n', [4,inf] ) ; 

fclose(fid) ;