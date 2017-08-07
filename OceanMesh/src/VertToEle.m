function [vtoe,nne]=VertToEle(t)
np = max(t(:)); ne = length(t); 
nne = zeros(np,1);                             
for ie = +1 : ne                                                            %--Go through once and find max nnz.
    for iv = +1 : +3
       vrtx = t(ie,iv); 
       nne(vrtx,1) = nne(vrtx,1) + 1; 
    end   
end
mnz  = max(nne);                                                           %--max number of non-zeros 
vtoe = zeros(np,mnz);                                                      %--vertex to element connectivity 
nne = zeros(np,1);                                                         %--number of neighboring elements 
for ie = +1 : ne
  for iv = +1 : +3 
     vrtx = t(ie,iv); 
     nne(vrtx,1) = nne(vrtx,1) +1; 
     vtoe(vrtx,nne(vrtx,1)) = ie; 
  end
end
nne  = nne'; 
vtoe = vtoe';

end