function [RefElObj] = InitRefOBJElVal( p )
%
% p is supposed to be a column vector
%
% Reference element:
%
%         REl.nNodes: 'number of nodes'
%       REl.RefNodes: 'nodal distribution'
%       RefVanderMat:  Van der Monde matrix V_{ij} = psi_{j}(x_{i})
%    InvRefVanderMat:  Inverse Van der Maonde matrix
%         RefMassMat:  Mass matrix <\phi_{i},phi_{j}>
%      InvRefMassMat:  Invese mass matrix
%              RefDr:  Derivative matrix wrt to r
%              RefDs:  Derivative matrix wrt to s
%              RefSr:  S_{r} - Weak derivative matrix wrt to r (weak form)
%              RefSs:  S_{s} - Weak derivative matrix wrt to s (weak form)
%                Nfp: 
%             nFaces: 
%               Edge: 
%       VecEdgeNodes: 
%            RefEmat: 
% 
%  See JSH, Nodal DG
%
%

RefElObj.p = p ; 
for i = 1: length(p)
    % 1. Compute reference nodes in [r, s] domain
    RefElObj.REl(i).nNodes = (p(i) + 1)*(p(i) + 2)/2 ; 
    [xx,yy] = Nodes2D( p(i) ) ;
    [r,s] = xytors( xx, yy ) ;
    RefElObj.REl(i).RefNodes = [r s] ;
    
    % 2. Compute a reference Vandermode 2D matrix
    RefElObj.REl(i).RefVanderMat = Vandermonde2D( p(i), r, s ) ;
    [nr,nc] = size(RefElObj.REl(i).RefVanderMat) ; 
    RefElObj.REl(i).InvRefVanderMat = RefElObj.REl(i).RefVanderMat\eye(nr,nc) ;
    
    % 3. Get a reference mass matrix
    MassMat = RefElObj.REl(i).RefVanderMat*(RefElObj.REl(i).RefVanderMat') ;
    [nr,nc] = size(MassMat) ; 
    
    % 3.1 A reference mass matrix
    RefElObj.REl(i).RefMassMat = MassMat\eye(nr,nc) ;
    % 3.2 Inverse of the reference mass matrix 
    RefElObj.REl(i).InvRefMassMat = MassMat ;  
    
    % 4. Derivative matrix
    [Dr,Ds] = Dmatrices2D( p(i), r, s, RefElObj.REl(i).RefVanderMat ) ;
    RefElObj.REl(i).RefDr = Dr ;
    RefElObj.REl(i).RefDs = Ds ; 
   
    % 4.1 Weak derivative matrix for weak form
    %     s_{i,j} = \int \phi_{j} \dfrac{\partial \phi_{i}}{\partial x}
    %     Note: for strong form, use S' for the weak derivative
    RefElObj.REl(i).RefSr = (RefElObj.REl(i).RefMassMat*Dr)' ; 
    RefElObj.REl(i).RefSs = (RefElObj.REl(i).RefMassMat*Ds)' ; 
    
    %
    % 4.2 Weak diffusion stiffness matrix
    %     k_{i,j} = \int \dfrac{\partial \phi_{i}}{\partial x} \dfrac{\partial \phi_{j}}{\partial x}
    %
    % Added on: Nov 21, 2014
    %
    RefElObj.REl(i).RefKrr = (Dr'*RefElObj.REl(i).RefMassMat*Dr) ; 
    RefElObj.REl(i).RefKss = (Ds'*RefElObj.REl(i).RefMassMat*Ds) ;
    %RefElObj.REl(i).RefKrs = (Dr'*RefElObj.REl(i).RefMassMat*Ds) ;
    
   %
   % 4.3. Second derivative matrices. Mainly y for SUPG.
   %
   [Drr,Dss,Drs] = D2matrices2D( p(i), r, s, RefElObj.REl(i).RefVanderMat ) ;
   RefElObj.REl(i).RefDrr = Drr ;
   RefElObj.REl(i).RefDss = Dss ;
   RefElObj.REl(i).RefDrs = Drs ;
    
    
    % 5. Compute necessary quantities assocaite with the reference faces
    %    of the reference lement
    nFaces = 3 ;
    Nfp = p(i) + 1 ;
    RefElObj.REl(i).Nfp = Nfp ; 
    RefElObj.REl(i).nFaces = nFaces ;
    
    % 5.1 Edge masks 
    fmask = zeros(Nfp,3) ;
    fmask(:,1) = find( abs(s + 1) < 1e-12 ) ; % face 1 
    fmask(:,2) = find( abs(r + s) < 1e-12 ) ; % face 2
    fmask(:,3) = find( abs(r + 1) < 1e-12 ) ; % face 3
    RefElObj.REl(i).Edge(1).Nodes = sort(fmask(:,1)) ;
    RefElObj.REl(i).Edge(2).Nodes = sort(fmask(:,2)) ;
    RefElObj.REl(i).Edge(3).Nodes = sort(fmask(:,3), 'descend') ;
    % RefElObj.REl(i).Edge(3).Nodes = sort(fmask(:,3)) ;
    
    % Edge nodes, counter clockwise 
    % RefElObj.REl(i).VecEdgeNodes = [fmask(:,1)' fmask(:,2)' fmask(:,3)']' ; 
    RefElObj.REl(i).VecEdgeNodes = [ RefElObj.REl(i).Edge(1).Nodes' ...
                                     RefElObj.REl(i).Edge(2).Nodes' ...
                                     RefElObj.REl(i).Edge(3).Nodes' ]' ; 
    
    % NOTE: this idea is taken from JSH                  
    % 6. compute mass matrix associcated with the edges     
    RefElObj.REl(i).RefEmat = zeros(length(xx), nFaces*Nfp ) ; 
    
    % Edge (1, 2, 3)
    FR = zeros(Nfp,3) ;
    FR(:,1) = r(RefElObj.REl(i).Edge(1).Nodes) ;
    FR(:,2) = r(RefElObj.REl(i).Edge(2).Nodes) ;
    FR(:,3) = s(RefElObj.REl(i).Edge(3).Nodes) ;
    for ie = 1: 3
        Vmat = Vandermonde1D( p(i), FR(:,ie) ) ;
        MassEdge = inv(Vmat*Vmat') ;
        % RefElObj.REl(i).RefEmat(fmask(:,ie),((ie - 1)*Nfp+1):(ie*Nfp)) = MassEdge ; 
        RefElObj.REl(i).RefEmat(RefElObj.REl(i).Edge(ie).Nodes,((ie - 1)*Nfp+1):(ie*Nfp)) = MassEdge ; 
    end
    %
    % Usage: Given Uel, fist set
    %   VecEl(RefElObj.REl(i).Edge(1).Nodes) = UEl(RefElObj.REl(i).Edge(1).Nodes)*GE.FDetJ(1)
    %   VecEl(RefElObj.REl(i).Edge(2).Nodes) = UEl(RefElObj.REl(i).Edge(2).Nodes)*GE.FDetJ(2)
    %   VecEl(RefElObj.REl(i).Edge(3).Nodes) = UEl(RefElObj.REl(i).Edge(3).Nodes)*GE.FDetJ(3)
    % and then one has 
    %   {\sum_{j} g_{j} \int \phi_{j} \phi_{i}} =
    %        RefElObj.REl.RefEmat*VecEl(RefElObj.REl(i).VecEdgeNodes)
    %-------------------
    % Note: the value of the leghth of each refence face is 2       
    %                                                                
    % Strickly, for face number 2, length(I_{2}) = 2*sqrt(2) but            
    % since I_{2} = sqrt(2) \bar{I}_{3}                              
    % and   I_{2}|_{phy} =  vol(I_{2}|{phy}) I_{2}/(2 sqrt(2))        
    %                    =  vol(I_{2}|{phy}) \bar{I_{2}}/2.
    % So that the value of the integral in the physical domain remains correct, if 
    % one uses length(I_{2}) = 2 and take
    %       I_{2} = \bar{I}_{2}
    %
end