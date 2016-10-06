function outname = writefort14_nodes( outfiname, VX, B)

% Read ADCIRC grid file fort.14

disp('read adcirc fort.14') ; 
tic
%if ( nargin == 0 ) 
%  finputname = 'fort.14'
%else
%  finputname = finame ;   
%end

if ( strcmpi(strtrim(outfiname),'fort.14') )
    disp('output file must not be "fort.14"') ; 
    return ; 
end

fid = fopen(outfiname,'w') ;
outname = outfiname ; 

title = 'interpolated bathymetry';
disp( title )  ;
fprintf(fid,'%s \n', title ) ; 

% print NP
NP = length(VX) ; 

% write NE, NP
fprintf(fid, '%d ! %s \n', NP, 'NP' ) ; 

% Improve write efficiency
% Vertices
pval = [ [1:1:NP]' VX B] ; 
fprintf( fid, '%10d %16.10f %16.10f %16.10f \n', pval' ) ; 

fclose(fid) ;

return

