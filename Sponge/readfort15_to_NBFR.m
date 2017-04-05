function f15dat = readfort15_to_NBFR( f15name )
%
%
% Readfort.15
%    
%

f15file = f15name ; 

fid = fopen( strtrim(f15file) ) ;

f15dat.rundes = fgetl(fid) ; % RUNDES
f15dat.runid = fgetl(fid) ;  % RUNID

% NFOVER
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.nfover = str2num(token) ;

% NABOUT
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.nabout = str2num(token) ;

% NSCREEN
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.nscreen = str2num(token) ;

% IHOT
line = fgetl(fid) ;
[token,res] = strtok(line) ;
f15dat.ihot = str2num(token) ;

% ICS
line = fgetl(fid) ;
[token,res] = strtok(line) ;
f15dat.ics = str2num(token) ;

% IM
f15dat.iden = [] ;
line = fgetl(fid) ; 
[token,res] = strtok(line) ;
f15dat.im = str2num(token) ;

% IDEN
if ( f15dat.im == 20 || f15dat.im == 30 ) 
    line = fgetl(fid) ;
    [token,res] = strtok(line) ; 
    
    f15dat.iden = str2num(token) ;  
end

% NOLIBF
token = readlinetoken( fid ) ;
f15dat.nolibf = str2num(token) ;

% NOLIFA
token = readlinetoken( fid ) ;
f15dat.nolifa = str2num(token) ;

% NOLICA
token = readlinetoken( fid ) ;
f15dat.nolica = str2num(token) ; 

% NOLICAT
token = readlinetoken( fid ) ; 
f15dat.nolicat = str2num( token ) ; 

% NWP
f15dat.nwp = 0 ;
token = readlinetoken( fid ) ;
f15dat.nwp = str2num( token ) ;
if ( f15dat.nwp > 0 ) 
    for l = 1: f15dat.nwp
        f15dat.AttrName(l).name = fgetl(fid) ; 
    end
end

% NCOR
token = readlinetoken( fid ) ;
f15dat.ncor = str2num( token ) ; 

% NTIP
token = readlinetoken( fid ) ; 
f15dat.ntip = str2num( token ) ;

% NWS
token = readlinetoken( fid ) ; 
f15dat.nws = str2num( token ) ;

% NRAMP
token = readlinetoken( fid ) ; 
f15dat.nramp = str2num( token ) ; 

% G
token = readlinetoken( fid ) ; 
f15dat.gravity = str2num( token ) ; 

% TAU0
token = readlinetoken( fid ) ;
f15dat.tau0 = str2num( token ) ; 

% Tau0FullDomainMin, Tau0FullDomainMax 
if ( abs(f15dat.tau0 + 5.D0) < 1e-10 )
    f15dat.tau0minmax =  readlinevec( fid ) ;
end

% DTDP
token = readlinetoken( fid ) ; 
f15dat.dtdp = str2num( token ) ; 

% STATIM
token = readlinetoken( fid ) ; 
f15dat.statim = str2num( token ) ; 

% REFTIM
token = readlinetoken( fid ) ; 
f15dat.reftim = str2num( token ) ; 

% WTIMINC
if ( ~mod(f15dat.nws - 1,10) || ... 
       ~mod(f15dat.nws - 11,10) ) 
  % token = readlinetoken( fid ) ; 
  % f15dat.wtimnc = str2num( token ) ; 
  f15dat.wtimnc = readlinevec( fid ) ;
end
  
% RNDY
token = readlinetoken( fid ) ; 
f15dat.rndy = str2num( token ) ;

% DRAMP
f15dat.dramp = readlinevec( fid ) ;

% A00, B00, C00
f15dat.a00b00c00 = readlinevec( fid ) ; 

% H0
f15dat.h0 = readlinevec( fid ) ;

% SLAM0, SFEA0
f15dat.slam = readlinevec( fid ) ; 

% TAU, CF, HBREAK, FTHETA, FGAMMA
if ( f15dat.nolibf <= 2 ) 
    f15dat.taucf = readlinevec( fid ) ;
end

% ESLM, ESLC
f15dat.elsm = [] ;
if ( f15dat.im <= 2 || f15dat.im == 10 )
    f15dat.elsm = readlinevec( fid ) ; 
end

% CORI
f15dat.cori = readlinescalar( fid ) ; 

% NTIF
f15dat.ntif = readlinescalar( fid ) ;

% Tidal potential
for k = 1: f15dat.ntif
    f15dat.tipotag(k).name = fgetl( fid ) ;
    f15dat.tipspec(k).val =  readlinevec( fid ) ; 
end
 
% NBFR
f15dat.nbfr = readlinescalar( fid ) ; 
for k = 1: f15dat.nbfr
    f15dat.bountag(k).name = fgetl( fid ) ; 
    f15dat.bounspec(k).val = readlinevec( fid ) ; 
end

fclose(fid) ;

    function vec = readlinevec( fid )
        msg = fgetl(fid) ; 
        
        idx = find(msg == ',') ; 
        if ( ~isempty(msg) ) 
            msg(idx) = ' ' ;
        end
        vec = sscanf(msg,'%f') ; 
    end

    function val = readlinescalar( fid )
        
       msg = fgetl(fid) ; 
       
       [ftoken,rem] = strtok(msg) ; 
       
       val = str2num(ftoken) ;
    end

    function token = readlinetoken( fid )
        
       msg = fgetl(fid) ; 
       [token,rem] = strtok(msg) ; 
        
    end
 
end