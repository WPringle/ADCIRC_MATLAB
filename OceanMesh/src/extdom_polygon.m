function [vs,ids,ide,error] = extdom_polygon( etbv, pv, iedbeg, ipsbeg, ipsend )
%
% Still very lousy
% 
ter = 1 ;
error = 0;

ipb  = ipsend ;
ipst = ipsbeg ;
idcur = iedbeg ;

sk = 1 ;
vs(sk,:) = pv(etbv(ipb,idcur),:) ;
ids(sk) = etbv(ipb,idcur) ;
while ( ter )
    if length(ipst) > 1
        error = 1;
        break;
    end
    ikk = etbv(ipst,idcur) ;
    
    ide(sk) = idcur ; 
    
    sk = sk + 1 ;
    vs(sk,:) = pv(ikk,:) ;
    ids(sk) = ikk ;
    
    itemp = find(~(etbv - ikk)) ; 
    iednext = ceil(itemp/2)  ;
   
    %
    if ~isempty(iednext(find((iednext - idcur))))
        idcur = iednext(find((iednext - idcur))) ; 
        ipst = find(etbv(:,idcur) - ikk) ;
    else
        idcur = iedbeg;
    end

    %
    
    if ( idcur == iedbeg ) 
        ter = 0 ;
    end
end