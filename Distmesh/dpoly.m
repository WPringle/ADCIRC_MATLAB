function [d,iloc]=dpoly(p,pv,ps)

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

% p are the mesh points
% pv is the bounded polygon
% ps is the points on polygon

% d is the distance from point, p to closest point on polygon, ps (d is
% negative if inside the bounded polygon, pv and positive if outside)
% iloc is the indice in polygon for closest point from p

np=size(p,1) ;
if ~isempty(ps)
    %nvs=size(pv,1)-1;
    nvs=size(ps,1)-1;
    memneed = np*nvs*8/10^9  ; % assume real number  
    if isempty(gcp('nocreate'))
        %disp('dsegment serial')
        maxmem = 4; % 4 Gb
        if ( memneed < maxmem ) % 8 Gb
           %ds=dsegment(p,pv) ;
           ds=dsegment(p,ps) ;
           [d,~]=min(ds,[],2);
        else
           csz = ceil(maxmem*10^9/(8*nvs)) ;

           ires = mod( np, csz) ;
           if ( ires )
               ichunck = [0:csz:np np] ;
           else
               ichunck = [0:csz:np] ; 
           end

           d = zeros(np,1) ;
           iloc = zeros(np,1) ; 
           nch = length(ichunck) ;
           for ic = 1: nch - 1
               % ic
               ibeg = ichunck(ic) + 1 ;
               iend = ichunck(ic+1)   ;

               pp = p(ibeg:iend,:) ; 

               %dsc = dsegment(pp,pv) ;
               dsc = dsegment(pp,ps) ;
               [d(ibeg:iend),iloc(ibeg:iend)] = min( dsc, [], 2) ;  
           end
        end
    else
        %disp('dsegment parallel')
        d = zeros(length(p),1); iloc = zeros(length(p),1);
        parfor ii = 1:np
            ds              = dsegment(p(ii,:),ps);
           [d(ii),iloc(ii)] = min(ds);
        end
    end
else
    % make distance very large as must be very far from any boundaries
    iloc = zeros(length(p),1);
    d = 1d8*ones(length(p),1);
end

if isempty(gcp('nocreate'))
    %disp('InPolygon serial')
    in = InPolygon(p(:,1),p(:,2),pv(:,1),pv(:,2));
else
    % Parallel calls to evaulate inpolygon function
    %disp('InPolygon parallel')
    in = zeros(length(p),1);
%     parfor ii = 1:np
%         in(ii) = InPolygon(p(ii,1),p(ii,2),pv(:,1),pv(:,2));
%         if mod(ii,1d5) == 0
%             disp(ii)
%         end
%     end
    Pool = gcp(); num_p = Pool.NumWorkers;
    for idx = 1:num_p
        ns = int64((idx-1)*np/num_p)+1;
        ne = int64(idx*np/num_p);
        f(idx) = parfeval(Pool,@InPolygon,1,...
                          p(ns:ne,1),p(ns:ne,2),pv(:,1),pv(:,2));
    end
    for idx = 1:num_p
        [idx_t, in_t] = fetchNext(f); % Get results into a cell array
        ns = int64((idx_t-1)*np/num_p)+1;
        ne = int64(idx_t*np/num_p);
        in(ns:ne) = in_t;
    end
end
d=(-1).^(in).*d;
%disp('dpoly out')
%

% MEXED

%function ds=donesegment(p,pv)
%
%e=ones(size(p,1),1);
%
%v=diff(pv,1);
%w=p-e*pv(1,:);
%
%c1=sum(w.*v(e,:),2);
%c2=sum(v(e,:).^2,2);
%
%ds=0*e;
%
%ix=c1<=0;
%ds(ix)=sqrt(sum((p(ix,:)-pv(1*ones(sum(ix),1),:)).^2,2));
%
%ix=c1>=c2;
%ds(ix)=sqrt(sum((p(ix,:)-pv(2*ones(sum(ix),1),:)).^2,2));
%
%ix=c1>0 & c2>c1;
%nix=sum(ix);
%if nix>0
%  Pb=ones(nix,1)*pv(1,:)+c1(ix)./c2(ix)*v;
%  ds(ix)=sqrt(sum((p(ix,:)-Pb).^2,2));
%end
