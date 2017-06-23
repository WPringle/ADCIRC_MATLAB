function [ampu,phiu,ampv,phiv,harminfo] = readfort54( finame )
%
%
fort54 = 'fort.54' ; 
if ( nargin >= 1 ) 
    fort54 = strtrim(finame) ;
end

fid = fopen(fort54) ;

% n frequency
harminfo.nfreq = fscanf(fid, '%f \n', 1 ) ; 

nfreq = harminfo.nfreq ; 
for i = 1: nfreq
    line = fgetl( fid )  ;
    harminfo.freqinfo(i).C = textscan(line,'%f %f %f %s')  ;
end

% line = fgetl( fid ) 
np = fscanf(fid, '%f \n', 1) ;


ampu = zeros(np,nfreq) ; 
phiu = zeros(np,nfreq) ; 

ampv = zeros(np,nfreq) ; 
phiv = zeros(np,nfreq) ; 

for i = 1: np
    id = fscanf(fid,'%f \n',1) ;
    for n = 1: nfreq
        val = fscanf(fid,'%f %f %f %f \n', [1 4]) ;
        ampu(i,n) = val(1) ;
        phiu(i,n) = val(2) ;
        ampv(i,n) = val(3) ;
        phiv(i,n) = val(4) ; 
    end
end

fclose(fid) ; 