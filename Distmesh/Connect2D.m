function [EToE, EToF] = Connect2D(EToV)
% function [EToE, EToF] = Connect2D(EToV)
% Purpose  : Build global connectivity arrays for grid based on
%            standard EToV input array from grid generator
% JSH

Nfaces = 3;

% Find number of elements and vertices
K = size(EToV,1); Nv = max(max(EToV));

% Create face to node connectivity matrix
TotalFaces = Nfaces*K;

% List of local face to local vertex connections
vn = [[1,2];[2,3];[1,3]];
% vn = [[1,2];[2,3];[3,1]] ; 

% tic
% % Build global face to node sparse array
% SpFToV = spalloc(TotalFaces, Nv, 2*TotalFaces);
% sk = 1;
% for k=1:K
%     for face=1:Nfaces
%         SpFToV( sk, EToV(k, vn(face,:))) = 1;
%         sk = sk+1;
%     end
% end
% toc

%tic
% 
%  Alternative version, much faster
%  Jan, 2016
vvn = reshape(vn',6,1) ;
icol = reshape(EToV(:,vvn)',6*K,1) ; 

irv = [1:3*K] ;
irv = reshape(irv,3,K)' ;

irow = zeros(K,6) ;

irow(:,1) = irv(:,1) ;
irow(:,3) = irv(:,2) ; 
irow(:,5) = irv(:,3) ; 
irow(:,2) = irv(:,1) ;
irow(:,4) = irv(:,2) ; 
irow(:,6) = irv(:,3) ; 
irow = reshape(irow',6*K,1) ;

val = icol*0 + 1 ; 
SpFToV = sparse(irow,icol,val) ; 

% max(max( SpFToV - SpFToV1 ))

%SpFToV(1,:)
%SpFToV1(1,:)
%
%toc


% Build global face to global face sparse array
SpFToF = SpFToV*SpFToV' - 2*speye(TotalFaces);

% Find complete face to face connections
[faces1, faces2] = find(SpFToF==2);

% Convert face global number to element and face numbers
element1 = floor( (faces1-1)/Nfaces )  + 1; face1 = mod( (faces1-1), Nfaces ) + 1;
element2 = floor( (faces2-1)/Nfaces )  + 1; face2 = mod( (faces2-1), Nfaces ) + 1;

% Rearrange into Nelements x Nfaces sized arrays
ind = sub2ind([K, Nfaces], element1, face1);

EToE = (1:K)'*ones(1,Nfaces); EToF = ones(K,1)*(1:Nfaces);
EToE(ind) = element2; EToF(ind) = face2 ;

return;
