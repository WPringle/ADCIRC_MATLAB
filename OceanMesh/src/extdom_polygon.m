function [poly,poly_idx,max_index,max_size] = extdom_polygon(bnde,pts,order)
% DESCRIPTION: Given a set of boundary edges of a singly- or multi- 
%              polygonal region, organize them in a winding order.
%
% INPUTS: 
%         bnde: the indices of each boundary edge as a nbnde x 2 matrix
%         pts:  the x,y locations of all the points in the region 
%               stored as an np x 2 matrix.
%         order:the order in which the traversal takes place
%               counter-clockwise (0) or clockwise (1). 
% OUTPUTS:
%          poly: the boundary of each enclosing polygon sorted in winding-order
%                poly is returned as a cell-array of length number of polys.
%      poly_idx: indices of the polygon coordinates in the same format as
%                poly
%     max_index: is the index into poly that is the largest
%      max_size: is the size of the largest poly
%
% kjr,UND,CHL,2017
%
%                                           TRAVERSAL METHOD
% Pick any unvisited edge segment [v_start,v_next] and add these vertices to the polygon loop.
% Find the unvisited edge segment [v_i,v_j] that has either v_i = v_next or v_j = v_next and add the other vertex (the one not equal to v_next) to the polygon loop.
% Reset v_next as this newly added vertex, mark the edge as visited and continue from 2.
% Traversal is done when we get back to v_start.
% NOTE: that the signed area will be positive if the vertices are
% oriented counterclockwise, and will be negative if it is oriented clockwise
bnde= unique(bnde,'rows');
active = true(size(bnde,1),1);
p = 0;
while any(active)
    p = p + 1;
    rn = find(active,1);
    v_start= bnde(rn,1);
    v_next = bnde(rn,2);
    active(rn) = false;
    temp  = pts(bnde(rn,:)',:);
    temp2 = bnde(rn,:)'; 
    k = 2;
    while v_next~=v_start
        [r,~]=find(v_next==bnde & active,1);
        tsel = bnde(r,:); sel=tsel(tsel~=v_next);
        k = k + 1;
        temp(k,:)= pts(sel,:);
        temp2(k,:)= sel;
        active(r) = false;
        v_next = sel;
    end
    poly{p}     = temp;
    poly_idx{p} = temp2;
    [area]=parea(poly{p}(:,1),poly{p}(:,2));
    if(order==0) % ccw
        if sign(area)<0
            poly{p} = flipud(poly{p});
            poly_idx{p} = flipud(poly_idx{p});
        end
    else % cw
        if sign(area)>0
            poly{p} = flipud(poly{p});
            poly_idx{p} = flipud(poly_idx{p});
        end
    end
end
[max_size, max_index] = max(cellfun('size', poly, 1));
end
% helper function, computes area of polygon
    function [area]=parea(x,y)
        n    = length(x);
        xp   = [x; x(1)];
        yp   = [y; y(1)];
        area = 0;
        for i = 1:n
            area = area + det([xp(i), xp(i+1); yp(i), yp(i+1)]);
        end
        area = 1/2*area;
    end