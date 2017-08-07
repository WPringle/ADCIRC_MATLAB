function [p,t] = distmesh2d(fd,fh,h0,bbox,p,pfix,itmax,plot_on,nscreen)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,P,PFIX,ITMAX,PLOT_ON,NSCREEN)
%
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      P:         Initial distribution of p if exist
%      PFIX:      Fixed node positions (NFIXx2)
%      PLOT_ON:   Flag to plot or not
%      NSCREEN:   Frequency to output temp mesh and info to screen

%   Modifications to Distmesh2d for quality check, removing low connectivity,
%   L0 calculation, Bessens-Heckbert force function, termination criterion,
%   clean-up.
%   William Pringle and Keith Roberts 2017
tic
it = +1 ; cl = 0;
geps=+.001*h0; deps=sqrt(eps)*h0;
ttol=+.1; dptol = +.001; Fscale=+1.2;
deltat = +0.1;
save_qual = zeros(itmax,2);
imp = 10; % number of iterations to do mesh improvements (delete/add)
if isempty(p)
    %% 1. Create initial distribution in bounding box (equilateral triangles)
    % Break it up into vertical blocks to avoid exceeding memory with ini.
    % point distribution
    noblks = +8;
    len = abs(bbox(1,1)-bbox(2,1));
    blklen = floor(len)/noblks;
    if len < +5; noblks = 1; end %--default back to one block if domain isn't too big.
    st = bbox(1,1) + blklen; ed = st + blklen;
    for blk = 1 : noblks
        if blk~=noblks
            [x,y]=meshgrid(st:h0:ed,...
                bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
            st = ed;
            ed = st + blklen;
        else
            [x,y]=meshgrid(st:h0:bbox(2,1),...
                bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
        end
        x(2:2:end,:)=x(2:2:end,:)+h0/+2;                % Shift even rows
        p1=[x(:),y(:)];                                 % List of node coordinates
        %% 2. Remove points outside the region, apply the rejection method
        p1=p1(feval(fd,p1,0)<geps,:);                   % Keep only d<0 points
        r0=1./feval(fh,p1).^+2;                         % Probability to keep point
        max_r0 = 1/h0^2;
        p1=p1(rand(size(p1,1),1) < r0/max_r0,:);        % Rejection method
        p = [p;p1];
    end
end
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end       
pfix=unique(pfix,'rows'); nfix=size(pfix,1);           % Remove duplicated nodes
p = [pfix; p];                                         % Prepend fix points
N = size(p,1);                                         % Number of points N
disp(['Number of initial points after rejection is ',num2str(N)]);
%% Iterate
pold = inf;                                             % For first iteration
if plot_on >= 1
    clf,view(2),axis equal,axis off
end
toc
fprintf(1,' ------------------------------------------------------->\n') ;
disp('Begin iterating...');
while 1
    tic
    if mod(it,nscreen) == 0
        disp(['Iteration =' num2str(it)]) ;
    end
    
    % 3. Retriangulation by the Delaunay algorithm
    if max(sqrt(sum((p-pold).^2,2))/h0)> ttol            % Any large movement?
        p = unique(p,'rows','stable');
        N = size(p,1); pold = p;                         % Save current positions
        t = delaunay_elim(p,fd,geps);                    % delaunay with elimination
        
        % 4. Describe each bar by a unique pair of nodes.
        bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
        bars=unique(sort(bars,2),'rows');                % Bars as node pairs
        
        % 5. Graphical output of the current mesh
        if plot_on >= +1 && (mod(it,nscreen)==0 || it == +1)
            cla,patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            title(['Iteration = ',num2str(it)]);
            drawnow
        end
    end
    
    % Getting element quality and check goodness
    tq = gettrimeshquan( p, t);
    mq_m = mean(tq.qm);
    mq_l = prctile(tq.qm,0.5);
    
    if  mod(it,imp) == imp-1 
        %[p,t] = fixmesh(p,t);
        %[p,t] = direct_smoother_lur(p,t); 
        %tqf  = Get_Finished_Quality(p,t);
        %minf = sort(tqf.qm);
        %minf = minf(20);
        %
        %minf = prctile(tqf.qm,0.05); %min(tqf.qm);
        %meanf = mean(tqf.qm);
        if mq_m > 0.95 && mq_l > 0.6
        %if minf > 0.5 && meanf > 0.975
              %[p,t] = remove_small_connectivity(p,t,fd,geps);
%             if ~isempty(nn)
%                 cl = cl + 1; % Only do this loop maximum of five times
%                 if cl < 5
%                     disp(['deleting ' num2str(length(nn)) ' due to small connectivity'])
%                     p(nn,:) = [];
%                     N = size(p,1); pold = inf;
%                     it = it + 1; 
%                     continue; 
%                 end
%            else
                disp('min finished quality of mesh is good enough, exit')
                close all;
                break;
%            end
%         else
%             N = size(p,1); pold = inf; it = it+1;
%             continue;
         end
    end
    
    % Saving a temp mesh
    if mod(it,nscreen) == 0
        disp(['Number of nodes is ' num2str(length(p))])
        disp(['Mean mesh quality is ' num2str(mq_m)])
        disp(['0.5 prctile mesh quality is ' num2str(mq_l)])
        save_qual(it,:) = [mq_m, mq_l];
        save('Temp_grid.mat','it','p','t');
    end
    
    % 6. Move mesh points based on bar lengths L and forces F
    barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
    hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2); % Ideal lengths
    L0=hbars*Fscale*nanmedian(L)/nanmedian(hbars);     % L0 = Desired lengths, our scale factor
    LN=L./L0;                                          % LN = Normalized bar lengths
    
    % Mesh improvements (deleting addition)
    if mod(it,imp)==0
        if save_qual(it-1,1) - save_qual(it-2,1)/deltat < +1e-2 
            % remove elements with small connectivity
            nn = get_small_connectivity(p,t);
            disp(['Deleting ' num2str(length(nn)) ' due to small connectivity'])
            % remove points that are too close.
            if any(L0 > +2*L)
                nn1 = setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix);
                disp(['Deleting ' num2str(length(nn1)) ' points too close together'])
                nn = unique([nn; nn1]);
            end
            % split too long edges
            il = find(LN>+1.5);
            if  ~isempty(il) 
                padd = (p(bars(il,1),:)+p(bars(il,2),:))/2;
                display(['Number of edges split : ',num2str(length(il))])
            end
            p(nn,:)= []; p = [padd; p];
            N = size(p,1); pold = inf; it = it+1;
            continue;
        end
    end
    
    F=(1-LN.^+4).*exp(-LN.^+4)./LN; % Bessens-Heckbert edge force
    Fvec =F*[1,1].*barvec;
    
    Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
    p=p+deltat*Ftot;                                   % Update node positions
    
    %7. Bring outside points back to the boundary
    d=feval(fd,p,0); ix=d>0; % Find points outside (d>0)
    dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],0)-d(ix))/deps; % Numerical
    dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],0)-d(ix))/deps; % gradient
    dgrad2 = dgradx.^+2 + dgrady.^+2;
    p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project
    
    % 8. Termination criterion: Exceed itmax
    it = it + 1 ;
    
    if ( it > itmax )
        nn = get_small_connectivity(p,t);
        if ~isempty(nn)
            cl = cl + 1; % Only do this loop maximum of five times
            if cl < 5
                disp(['deleting ' num2str(length(nn)) ' due to small connectivity'])
                p(nn,:) = [];
                N = size(p,1); pold = inf;
                continue;
            end
        end
        disp('too many iterations, exit')
        close all;
        break ;
    end
    
    % or all interior nodes move less than dptol (scaled)
    if ( max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0) < dptol )
        nn = get_small_connectivity(p,t);
        if ~isempty(nn)
            cl = cl + 1; % Only do this loop maximum of five times
            if cl < 5
                disp(['deleting ' num2str(length(nn)) ' due to small connectivity'])
                p(nn,:) = [];
                N = size(p,1); pold = inf;
                continue;
            end
        end
        disp('no movement of nodes, exit')
        close all;
        break;
    end
    toc
end
disp('Finished iterating...');
fprintf(1,' ------------------------------------------------------->\n') ;
%% Clean up
[p,t] = fixmesh(p,t);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary subfunctions %
%%%%%%%%%%%%%%%%%%%%%%%%%%

    function t = delaunay_elim(p,fd,geps)
        p_s = p - mean(p);
        t = DelaunayTri(p_s);
        %t = delaunayn(p_s);
        pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/+3;  % Compute centroids
        t = t(feval(fd,pmid,0)<-geps,:);                  % Keep interior triangles
    end

    function nn = get_small_connectivity(p,t)
        % Get node connectivity (look for 4)
        [~, nn] = VertToEle(t);
        %nn = NodeConnect2D( t );
        % Make sure they are not boundary nodes
        bdbars = extdom_edges2(t, p);
        bdnodes = unique(bdbars(:));
        I = find(nn <= 4);
        nn = setdiff(I',bdnodes);
        return;
    end

    function [p,t] = remove_small_connectivity(p,t,fd,geps)
        for iter = 1:5
            % Get node connectivity (look for 4)
            nn2 = NodeConnect2D( t );
            % Make sure they are not boundary nodes
            bdbars = extdom_edges(t, p);
            bdnodes = unique(bdbars(:));
            I = find(nn2 <= 4);
            nn2 = setdiff(I',bdnodes);
            if isempty(nn2); break; end
            p(nn2,:) = []; 
            t = delaunay_elim(p,fd,geps);
        end
        return; 
    end    

    function tq = Get_Finished_Quality(p,t)
        % Before exiting we fix bad edges and smooth mesh
        [p,t] = Fix_bad_edges_and_mesh(p,t,0);
        %[p,t] = smoothmesh(p,t,1000); 
        [p,t] = direct_smoother_lur(p,t); 
        tq = gettrimeshquan( p, t);
        minq = sort(tq.qm);
        tria = tarea(p,t);
        pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/+3;  % Compute centroids
    end

    function [ar] = tarea(p,t)
        %Compute triangle area
        %
        it1=t(:,1); it2=t(:,2); it3=t(:,3);
        x21=p(it2,1)-p(it1,1); y21=p(it2,2)-p(it1,2);
        x31=p(it3,1)-p(it1,1); y31=p(it3,2)-p(it1,2);
        x32=p(it3,1)-p(it2,1); y32=p(it3,2)-p(it2,2);
        ar=(x21.*y31-y21.*x31)/2;
    end

    function[ncon,bn]=ng(t,np)
        % Returns the nodal graph (ncon) and the idicies of the
        % nodes on the boundary.
        ncon = min(sparse(t(:,1),t(:,2),1,np,np)+sparse(t(:,2),t(:,3),1,np,np)+sparse(t(:,3),t(:,1),1,np,np),+1);
        ncon = min(ncon+ncon',+1);
        % this finds the boundary points from the nodal graph, whose indexes are stored in bn
        B = ncon^2.*ncon==+1;
        bn = find(sum(B,2)>0);
    end

    function[bnei]=getBouNei(bn,ncon)
        % Returns the adj boundary neighbors to boundary node bn
        n=0*bn; bnei=zeros(length(bn),2);
        for i = +1 : length(bn)
            n=find(ncon(bn(i),:)==1);
            bnei(i,1:2)=intersect(n,bn);
        end
    end

end