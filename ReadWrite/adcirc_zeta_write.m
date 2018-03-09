function adcirc_zeta_write(t,x,y,e,d,z,file)
% adcirc_zeta_write  Write an adcirc .63.nc file for GMT/FigureGen
%
% Uses built-in NetCDF capability (MATLAB R2008b or later) to 
% write a COARDS-compliant netCDF grid file
% Duplicates (some) functionality of the program grdwrite (which requires
% compilation as a mexfile-based function on each architecture) using
% Matlab 2008b (and later) built-in NetCDF functionality
% instead of GMT libraries.
%
% adcirc_zeta_write(t,x,y,e,d,z,file) will create a grid file containing the
% data in the matrix Z.  X and Y should be either vectors with
% dimensions that match the size of Z or two-component vectors
% containing the max and min values for each.
%
% See also GRDREAD2, GRDINFO2

% For more information on GMT grid file formats, see:
% http://www.soest.hawaii.edu/gmt/gmt/doc/gmt/html/GMT_Docs/node70.html
% Details on Matlab's native netCDF capabilities are at:
% http://www.mathworks.com/access/helpdesk/help/techdoc/ref/netcdf.html

% GMT (Generic Mapping Tools, <http://gmt.soest.hawaii.edu>)
% was developed by Paul Wessel and Walter H. F. Smith

% Kelsey Jordahl
% Marymount Manhattan College
% http://marymount.mmm.edu/faculty/kjordahl/software.html

% Time-stamp: <Tue Jul 19 16:28:24 EDT 2011>
 
% Version 1.1.2, 19-Jul-2011
% Available at MATLAB Central
% <http://www.mathworks.com/matlabcentral/fileexchange/26290-grdwrite2>

if nargin < 4,
  help(mfilename);
  return,
end

% check for appropriate Matlab version (>=7.7)
V=regexp(version,'[ \.]','split');
if (str2num(V{1})<7) || (str2num(V{1})==7 && str2num(V{2})<7),
  ver
  error('grdread2: Requires Matlab R2008b or later!');
end

ncid = netcdf.create(file, 'NC_SHARE');
if isempty(ncid),
  return,
end

% set descriptive variables
conv='COARDS/CF-1.0';
title=file;
history='File written by MATLAB function adcirc_zeta_write.m';
desc=['Created ' datestr(now)];
vers='49';                            

% check X and Y
if (~isvector(x) || ~isvector(y)),
  error('X and Y must be vectors!');
end
if (length(x) ~= size(z,1)),    % number of columns don't match size of x
  minx=min(x); maxx=max(x);
  dx=(maxx-minx)/(size(z,3)-1);
  x=minx:dx:maxx;                       % write as a vector
end
if (length(y) ~= size(z,1)),    % number of columns don't match size of y
  miny=min(y); maxy=max(y);
  dy=(maxy-miny)/(size(z,2)-1);
  y=miny:dy:maxy;                       % write as a vector
end    

% match Matlab class to NetCDF data type
switch class(z)
 case 'single'
  nctype='NC_FLOAT';
  nanfill=single(NaN);
 case 'double'
  nctype='NC_DOUBLE';
  nanfill=double(NaN);
 case 'int8'
  nctype='NC_BYTE';
  nanfill=intmin(class(z));
  disp(['Warning: ''No data'' fill value set to ' num2str(nanfill)])
 case 'int16'
  nctype='NC_SHORT';
  nanfill=intmin(class(z));
  disp(['Warning: ''No data'' fill value set to ' num2str(nanfill)])
 case 'int32'
  nctype='NC_INT';
  nanfill=intmin(class(z));
  disp(['Warning: ''No data'' fill value set to ' num2str(nanfill)])
 otherwise
  error(['Don''t know how to handle data of class ''' class(z) '''.  Try converting to a supported data type (int8, int16, int32, single or double).'])
end

% global
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions',conv);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'title',title);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'history',history);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'description',desc);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'FigGen_version',vers);

zunits = 'metres';
zname = 'free surface level above geoid';

% T
dimid = netcdf.defDim(ncid,'time',length(t));
varid = netcdf.defVar(ncid,'time','double',dimid);
netcdf.putAtt(ncid,varid,'long_name','model time');
netcdf.putAtt(ncid,varid,'actual_range',[min(t) max(t)]);
netcdf.putAtt(ncid,varid,'units','seconds since 2011-03-11 14:46 +9:00');
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,t); 
% X
netcdf.reDef(ncid);
dimid = netcdf.defDim(ncid,'node',length(x));
varid = netcdf.defVar(ncid,'x','double',dimid);
netcdf.putAtt(ncid,varid,'long_name','Easting');
netcdf.putAtt(ncid,varid,'actual_range',[min(x) max(x)]);
netcdf.putAtt(ncid,varid,'units','metres in UTM coordinates using old Tokyo datum');
netcdf.putAtt(ncid,varid,'positive','east');
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,x);
% Y
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'y','double',dimid);
netcdf.putAtt(ncid,varid,'long_name','Northing');
netcdf.putAtt(ncid,varid,'actual_range',[min(y) max(y)]);
netcdf.putAtt(ncid,varid,'units','metres in UTM coordinates using old Tokyo datum');
netcdf.putAtt(ncid,varid,'positive','north');
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,y);
% ELEMENT
netcdf.reDef(ncid);
netcdf.defDim(ncid,'nele',length(e));
netcdf.defDim(ncid,'nvertex',3);
varid = netcdf.defVar(ncid,'element','int',[3 2]);
netcdf.putAtt(ncid,varid,'long_name','element');
netcdf.putAtt(ncid,varid,'standard_name','face_node_connectivity');
netcdf.putAtt(ncid,varid,'start_index',1);
netcdf.putAtt(ncid,varid,'units','nondimensional');
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,e);
% D
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'depth','double',1);
netcdf.putAtt(ncid,varid,'long_name','depth below geoid');
netcdf.putAtt(ncid,varid,'actual_range',[min(d) max(d)]);
netcdf.putAtt(ncid,varid,'coordinates','y x');
netcdf.putAtt(ncid,varid,'location','node');
netcdf.putAtt(ncid,varid,'units',zunits);
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,d);
% Z
netcdf.reDef(ncid);
varid = netcdf.defVar(ncid,'zeta',nctype,[1 0]);
netcdf.putAtt(ncid,varid,'long_name',zname);
netcdf.putAtt(ncid,varid,'actual_range',[min(z(:)) max(z(:))]);
netcdf.putAtt(ncid,varid,'coordinates','time y x');
netcdf.putAtt(ncid,varid,'location','node');
netcdf.putAtt(ncid,varid,'units',zunits);
netcdf.putAtt(ncid,varid,'_FillValue',-99999);%nanfill
netcdf.endDef(ncid);
netcdf.putVar(ncid,varid,z);
% close file
netcdf.close(ncid);
