function [pm,tm]=Grid_Merge(mesh1,mesh2)
% Merge 2-d simplex meshes compromised of np vertices and nt triangles
% contained in the mat files  mesh1 with mesh2
% Uses a simplified form of the algorithm described in the pulbication
% Connectivity Oblivious Merging of Triangulations

% INPUTS:
% mesh1: filename of base mesh.
% mesh2: filename of mesh to be merged into the base mesh.

% OUPUTS:
% pm: np x 2 coordinates of vertices of the merged mesh.
% tm: nt x 3 matrix representing the 2-d simplices.
% kjr, und, chl, sept. 2017 Version 1.0.
mesh{1}=load(mesh1);
p1=mesh{1}.p;
t1=mesh{1}.t; 

mesh{2}=load(mesh2);
p2=mesh{2}.p;
t2=mesh{2}.t; 

% merge two meshes p1 replaces the overlap in p2.
disp('Merging meshes...')
[pm,tm]=merge_meshes(p1,p2);

disp('Pruning outer triangles...')
% prune triangles outside both domains. 
pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/+3;

% form outer polygon 1 
bnde=extdom_edges2(t1,p1); 
poly1=extdom_polygon(bnde,p1,1); 
k=0; poly_vec1=[];
for i = 1 : length(poly1) 
   for ii = 1 : length(poly1{i}) 
      k = k + 1; 
      poly_vec1(k,:) = poly1{i}(ii,:);
   end
   k = k + 1; 
   poly_vec1(k,:) = [NaN,NaN];
end
[edges]=Get_poly_edges(poly_vec1);
in1=inpoly(pmid,poly_vec1,edges);

% form outer polygon 2
bnde=extdom_edges2(t2,p2); 
poly2=extdom_polygon(bnde,p2,1); 
k=0; poly_vec2=[];
for i = 1 : length(poly2) 
   for ii = 1 : length(poly2{i}) 
      k = k + 1; 
      poly_vec2(k,:) = poly2{i}(ii,:);
   end
   k = k + 1; 
   poly_vec2(k,:) = [NaN,NaN];
end
[edges]=Get_poly_edges(poly_vec2);
in2=inpoly(pmid,poly_vec2,edges);

tm(~in1 & ~in2,:) = [];

disp('Cleaning up')
% remove triangles w/ small connectivity (valency < 4).
% and then reproject onto the plane
nn = get_small_connectivity(pm,tm); pm(nn,:)= [];
[pm,tm]=remove_overlap(pm);

% prune triangles outside the domain agaiin.
pmid = (pm(tm(:,1),:)+pm(tm(:,2),:)+pm(tm(:,3),:))/+3;

% form outer polygon 1
[edges]=Get_poly_edges(poly_vec1);
in1=inpoly(pmid,poly_vec1,edges);

% form outer polygon 2
[edges]=Get_poly_edges(poly_vec2);
in2=inpoly(pmid,poly_vec2,edges);

tm(~in1 & ~in2,:) = [];


% clean up some more to avoid non-unique boundary edges.
[pm,tm]=Fix_bad_edges_and_mesh(pm,tm,0);
tq=gettrimeshquan(pm,tm);
tm(abs(tq.qm)<0.10,:)=[];

disp('Smoothing mesh...')
% Laplacian smoother to improve quality of merger region.
[pm,tm]=smoothmesh(pm,tm);

figure; simpplot(pm,tm);
title('merged');

return

    function nn = get_small_connectivity(p,t)
        % Get node connectivity (look for 4)
        [~, nn] = VertToEle(t);
        % Make sure they are not boundary nodes
        bdbars = extdom_edges2(t, p);
        bdnodes = unique(bdbars(:));
        I = find(nn <= 4);
        nn = setdiff(I',bdnodes);
        return;
    end

    function [p,t]=remove_overlap(p)
        % remove overlapping elements by projecting onto the 
        % lower convex hull and reprojecting to the plane.
        p(:,3) = p(:,1).^2 + p(:,2).^2;
        t=convhull(p(:,1),p(:,2),p(:,3));
        nt=length(t);
        r0 = p(t(:,1),:);
        r1 = p(t(:,2),:);
        r2 = p(t(:,3),:);
        
        vec1 = r1 - r0;
        vec2 = r2 - r0;
        
        % compute the face normal
        n = cross(vec1,vec2);
        e = [0,0,1]; % unit normal of plane
        islower = dot(n,repmat(e,[nt,1]),2) < eps;
        
        t(~islower,:)=[];
        p(:,3)=[];
    end

    function [pm,tm]=merge_meshes(p1,p2)
        % merge triangulation 1 with triangulation2 2
        % the vertices of triangulation 1 are in p1 (no. vert x 2)
        % the vertices of triangulation 2 are in p2 (no. vert x 2)
        % 1. project onto paraboloid
        % 2. determine lower convex hull
        % 3. project lower hull onto plane
        
        z = p1(:,1).^2 + p1(:,2).^2 + 10e-4; % weight is constant for DT
        z2= p2(:,1).^2 + p2(:,2).^2 + 9.0e-4;  % weight for t2 so it appears on hull
        
        xm = [p1(:,1);p2(:,1)]; ym = [p1(:,2);p2(:,2)]; zm = [z;z2];
        pm = [xm,ym,zm];
        tm=convhull(xm,ym,zm);
        nt=length(tm);
        
        % r0 is anchor point
        r0 = pm(tm(:,1),:);
        r1 = pm(tm(:,2),:);
        r2 = pm(tm(:,3),:);
        
        vec1 = r1 - r0;
        vec2 = r2 - r0;
        
        % compute the face normal
        n = cross(vec1,vec2);
        e = [0,0,1]; % unit normal of plane
        islower = dot(n,repmat(e,[nt,1]),2) < eps;
        
        tm(~islower,:)=[];
        pm(:,3)=[];
    end
end
