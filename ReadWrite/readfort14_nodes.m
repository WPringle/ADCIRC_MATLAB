function [VX, B] = readfort14_nodes( finame ) 
% Read ADCIRC grid file fort.14 nodes only

tic
if ( nargin == 0 ) 
  finputname = 'fort.14'
else
  finputname = finame ;   
end

fid = fopen(finputname) ;

agrid = fgetl(fid) ;
disp(agrid) ;
%title = agrid ; 

N = fscanf(fid,'%g %g',2) ;
%
% Nov 15, 2012, improve reading efficiency
Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;
%
% Arrange it to a Nodal DG input
VX = Val(:,2:3) ;
B  = Val(:,4) ;

fclose(fid) ;

return

