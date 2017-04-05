function fid = writefort15_to_BND( f15out, f15dat ) 
%
if ( strcmp(strtrim(f15out), 'fort.15') )
    disp('Error: an output file name must not be fort.15') ;
    return ;
end

fid = fopen( f15out, 'w' ) ;

fprintf(fid, '%s\n', f15dat.rundes ) ; % RUNDES
fprintf(fid, '%s\n', f15dat.runid  ) ;  % RUNID

% NFOVER
% line = fgetl(fid) ; 
% [token,res] = strtok(line) ;
fprintf(fid, '%d        \t ! NFOVER \n', f15dat.nfover ) ;

% NABOUT
% line = fgetl(fid) ; 
% [token,res] = strtok(line) ;
fprintf(fid, '%d        \t ! NABOUT \n', f15dat.nabout ) ;

% NSCREEN
% line = fgetl(fid) ; 
% [token,res] = strtok(line) ;
fprintf(fid, '%d        \t ! NSCREEN \n', f15dat.nscreen ) ;

% IHOT
% line = fgetl(fid) ;
% [token,res] = strtok(line) ;
fprintf(fid, '%d        \t ! IHOT \n', f15dat.ihot ) ;

% ICS
% line = fgetl(fid) ;
% [token,res] = strtok(line) ;
fprintf(fid, '%d         \t ! ICS \n', f15dat.ics ) ;

% IM
%f15dat.iden = [] ;
%line = fgetl(fid) ; 
%[token,res] = strtok(line) ;
fprintf(fid,'%d          \t ! IM \n', f15dat.im ) ;

% IDEN
if ( f15dat.im == 20 || f15dat.im == 30 ) 
    % line = fgetl(fid) ;
    % [token,res] = strtok(line) ; 
    % 
    % f15dat.iden = str2num(token) ;
    fprintf( fid, '%d           \t ! IDEN \n', f15dat.iden ) ; 
end

% NOLIBF
% token = readlinetoken( fid ) ;
fprintf( fid, '%d            \t ! NOLIBF \n', f15dat.nolibf ) ;

% NOLIFA
% token = readlinetoken( fid ) ;
fprintf( fid, '%d            \t ! NOLIFA \n', f15dat.nolifa ) ;

% NOLICA
% token = readlinetoken( fid ) ;
fprintf( fid, '%d           \t ! NOLICA \n', f15dat.nolica ) ; 

% NOLICAT
% token = readlinetoken( fid ) ; 
fprintf( fid, '%d           \t ! NOLICAT \n', f15dat.nolicat ) ; 

% NWP
% f15dat.nwp = 0 ;
% token = readlinetoken( fid ) ;

 fprintf(fid, '%d           \t ! NWP \n', f15dat.nwp ) ;
if ( f15dat.nwp > 0 ) 
    f15dat.nwp
    for l = 1: f15dat.nwp
        fprintf(fid, '%s\n', f15dat.AttrName(l).name ) ; 
    end
end

% NCOR
% token = readlinetoken( fid ) ;
fprintf(fid, '%d   \t \t ! NCOR \n', f15dat.ncor ) ; 

% NTIP
% token = readlinetoken( fid ) ; 
fprintf(fid, '%d   \t \t ! NTIP \n', f15dat.ntip ) ;

% NWS
% token = readlinetoken( fid ) ; 
fprintf(fid, '%d   \t \t ! NWS \n', f15dat.nws ) ;

% NRAMP
% token = readlinetoken( fid ) ; 
fprintf(fid, '%d    \t \t ! NRAMP \n', f15dat.nramp ) ;

% G
% token = readlinetoken( fid ) ; 
fprintf(fid, '%e    \t ! G  \n', f15dat.gravity ) ; 

% TAU0
% token = readlinetoken( fid )
fprintf(fid, '%16.9e   \t \t ! TAU0 \n', f15dat.tau0 ) ; 

% Tau0FullDomainMin, Tau0FullDomainMax 
if ( abs(f15dat.tau0 + 5.D0) < 1e-10 )
    fprintf(fid, '%e %e   \t ! Tau0FullDomainMin, Tau0FullDomainMax \n', f15dat.tau0minmax ) ;
end

% DTDP
% token = readlinetoken( fid ) ; 
fprintf( fid, '%16.9e    \t \t ! DTDP \n', f15dat.dtdp ) ; 

% STATIM
% token = readlinetoken( fid ) ; 
fprintf( fid, '%e  \t \t ! STATIM \n',  f15dat.statim ) ; 

% REFTIM
% token = readlinetoken( fid ) ; 
fprintf( fid, '%e   \t \t ! REFTIM \n', f15dat.reftim ) ; 

% WTIMINC
if ( ~mod(f15dat.nws - 1,10) || ... 
       ~mod(f15dat.nws - 11,10) ) 
  % token = readlinetoken( fid ) ; 
  % f15dat.wtimnc = str2num( token ) ; 
  fprintf( fid, '%e  ', f15dat.wtimnc ) ;
  fprintf( fid, '    \t \t ! WTMINC \n' ) ; 
end
  
% RNDY
% token = readlinetoken( fid ) ; 
fprintf( fid, '%g   \t \t ! RNDY \n', f15dat.rndy ) ;

% DRAMP
fprintf( fid, '%d   \t \t ! DRAMP \n',  f15dat.dramp ) ;

% A00, B00, C00
fprintf( fid, '%16.9e %16.9e %16.9e   \t  ! A00, B00, C00 \n', f15dat.a00b00c00 ) ; 

% H0
len = length(f15dat.h0) ;
for k = 1: len
  fprintf( fid, '%16.9e ', f15dat.h0(k) ) ;
end
fprintf( fid, '   \t ! H0 \n' ) ; 

% SLAM0, SFEA0
fprintf( fid, '%16.9e %16.9e  \t \t ! SLAM0, SFEA0 \n', f15dat.slam ) ; 

% TAU, CF, HBREAK, FTHETA, FGAMMA
if ( f15dat.nolibf <= 2 ) 
    fprintf( fid, '%16.9e ', f15dat.taucf ) ;
    fprintf( fid, '    \t ! TAU \n' ) ;
end

% ESLM, ESLC
if ( f15dat.im <= 2 || f15dat.im == 10 )
    fprintf( fid, '%16.9e ', f15dat.elsm ) ; 
    fprintf( fid, '    \t \t ! ELSM \n' ) ; 
end

% CORI
fprintf( fid, '%16.9e   \t ! CORI \n', f15dat.cori ) ; 

% NTIF
fprintf( fid, '%d   \t \t ! NTIF \n', f15dat.ntif ) ;

% Tidal potential
for k = 1: f15dat.ntif
    fprintf( fid, '%s \n',  f15dat.tipotag(k).name ) ;
    fprintf( fid, '%15.8e ', f15dat.tipspec(k).val ) ;
    fprintf( fid, '   \t !  TPK, AMIGT, ETRF, FFT, FACET \n' ) ; 
end
 
% NBFR
fprintf( fid, '%d   \t \t ! NBFR \n', f15dat.nbfr ) ; 
for k = 1: f15dat.nbfr
    fprintf( fid, '%s  \n', f15dat.bountag(k).name ) ; 
    fprintf( fid, '%15.8e %15.8e %15.8e \n', f15dat.bounspec(k).val ) ; 
end

% Open boundary harmonic forcing  
for k = 1: f15dat.nbfr
    fprintf(fid, '%s \n', f15dat.bountag(k).name  ) ; 
    
    % val = fscanf(fid, '%f %f \n' ) ; 
    fprintf(fid, '%16.9e %16.9e \n', f15dat.opeemoefa(k).val' ) ; 
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