function [amp,phi,harminfo,coors] = readfort53( finame )
%
%
fort53 = 'fort.53' ; 
if ( nargin >= 1 ) 
    fort53 = strtrim(finame) ;
end

fid = fopen(fort53) ;

% n frequency
harminfo.nfreq = fscanf(fid, '%f \n', 1 ) ; 

nfreq = harminfo.nfreq ;
for i = 1: nfreq
    line = fgetl( fid )  ;
    harminfo.freqinfo(i).C = textscan(line,'%f %f %f %s')  ;
end

% line = fgetl( fid ) 
np = fscanf(fid, '%f \n', 1) ;

amp = zeros(np,nfreq) ; 
phi = zeros(np,nfreq) ; 
coors = zeros(np,2) ; 

for i = 1: np
    line = fgetl( fid ) ; 
    id = sscanf(line,'%f') ;
    if ( length(id) > 1 )
        coors(i,:) = id(2:3) ; 
    else
        coors(i,1) = id ; 
    end
    for n = 1: nfreq
        val = fscanf(fid,'%f %f \n', [1 2]) ;
        amp(i,n) = val(1) ;
        phi(i,n) = val(2) ;
    end
end

fclose(fid) ; 