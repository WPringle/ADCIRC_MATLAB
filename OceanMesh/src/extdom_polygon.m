function [poly,max_index,max_size] = extdom_polygon(bnde,pts,order)
% DESCRIPTION: Given a set of boundary edges from a 2-d triangulation of 
%              points, organize them in a winding order.
%
% INPUTS: 
%         bnde: the vertices of each boundary edge as a nbnde x 2 matrix
%         pts:  the x,y locations of all the points in the triangulation 
%               stored as an np x 2 matrix.
%         order:the order in which the traversal takes place
%               counter-clockwise (0) or clockwise (1). 
% OUTPUTS:
%          poly: the boundary of each enclosing polygon sorted in winding-order
%                poly is returned as a cell-array of length number of polys.
%          max_index: is the index into poly that is the largest
%          max_size: is the size of the largest poly
%
% kjr,UND,CHL,2017
%
%                                           TRAVERSAL METHOD
% Pick any unvisited edge segment [v_start,v_next] and add these vertices to the polygon loop.
% Find the unvisited edge segment [v_i,v_j] that has either v_i = v_next or v_j = v_next and add the other vertex (the one not equal to v_next) to the polygon loop.
% Reset v_next as this newly added vertex, mark the edge as visited and continue from 2.
% Traversal is done when we get back to v_start.
% NOTE: that the signed area will be positive if the vertices V0V1V2 are
% oriented counterclockwise, and will be negative if the triangle is oriented clockwise
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
    k = 2;
    while v_next~=v_start
        [r,c]=find(v_next==bnde & active,1);
        tsel = bnde(r,:); sel=tsel(tsel~=v_next);
        k = k + 1;
        temp(k,:)= pts(sel,:);
        active(r) = false;
        v_next = sel;
    end
    poly{p} = temp;
    [area]=parea(poly{p}(:,1),poly{p}(:,2));
    if(order==0) % ccw
        if sign(area)<0
            poly{p} = flipud(poly{p});
        end
    else % cw
        if sign(area)>0
            poly{p} = flipud(poly{p});
        end
    end
end
[max_size, max_index] = max(cellfun('size', poly, 1));
end
% helper function
    function [area]=parea(x,y)
        n = length(x);
        xp = [x; x(1)];
        yp = [y; y(1)];
        area = 0;
        for i = 1:n
            area = area + det([xp(i), xp(i+1); yp(i), yp(i+1)]);
        end
        area = 1/2*area;
    end