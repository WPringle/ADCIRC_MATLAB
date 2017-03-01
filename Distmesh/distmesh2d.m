function [p,t]=distmesh2d(fd,fh,h0,bbox,p,pfix,itmax,plot_on)
%DISTMESH2D 2-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESH2D(FD,FH,H0,BBOX,PFIX,plot_on)
%
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      FCI:       Change interpolant function for points with small
%                 interiori angle
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      P:         Initial distribution of p if exist
%      PFIX:      Fixed node positions (NFIXx2)
%      ITMAX      Maximum number of iterations allowed
%      PLOT_ON:   Flag to plot or not
%
%   Example: (Uniform Mesh on Unit Circle)
%      fd=@(p) sqrt(sum(p.^2,2))-1;
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-1,-1;1,1],[]);
%
%   Example: (Rectangle with circular hole, refined at circle boundary)
%      fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
%      fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
%      [p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
%
%   Example: (Polygon)
%      pv=[-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
%          1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5];
%      [p,t]=distmesh2d(@dpoly,@huniform,0.1,[-1,-1; 2,1],pv,pv);
%
%   Example: (Ellipse)
%      fd=@(p) p(:,1).^2/2^2+p(:,2).^2/1^2-1;
%      [p,t]=distmesh2d(fd,@huniform,0.2,[-2,-1;2,1],[]);
%
%   Example: (Square, with size function point and line sources)
%      fd=@(p) drectangle(p,0,1,0,1);
%      fh=@(p) min(min(0.01+0.3*abs(dcircle(p,0,0,0)), ...
%                   0.025+0.3*abs(dpoly(p,[0.3,0.7; 0.7,0.5]))),0.15);
%      [p,t]=distmesh2d(fd,fh,0.01,[0,0;1,1],[0,0;1,0;0,1;1,1]);
%
%   Example: (NACA0012 airfoil)
%      hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4;
%      a=.12/.2*[0.2969,-0.1260,-0.3516,0.2843,-0.1036];
%
%      fd=@(p) ddiff(dcircle(p,circx,0,circr),(abs(p(:,2))-polyval([a(5:-1:2),0],p(:,1))).^2-a(1)^2*p(:,1));
%      fh=@(p) min(min(hlead+0.3*dcircle(p,0,0,0),htrail+0.3*dcircle(p,1,0,0)),hmax);
%
%      fixx=1-htrail*cumsum(1.3.^(0:4)');
%      fixy=a(1)*sqrt(fixx)+polyval([a(5:-1:2),0],fixx);
%      fix=[[circx+[-1,1,0,0]*circr; 0,0,circr*[-1,1]]'; 0,0; 1,0; fixx,fixy; fixx,-fixy];
%      box=[circx-circr,-circr; circx+circr,circr];
%      h0=min([hlead,htrail,hmax]);
%
%      [p,t]=distmesh2d(fd,fh,h0,box,fix);

%
%   See also: MESHDEMO2D, DISTMESHND, DELAUNAYN, TRIMESH.

%   distmesh2d.m v1.1
%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.
tic
it = 1 ;
% dxx = (bbox(2,1) - bbox(1,1)) ;
% dyy = (bbox(2,2) - bbox(1,2)) ;
%
% deltat = h0/20 
% itmax = max(dxx/deltat, dyy/deltat) ; 

dptol=.0008; ttol=.1; Fscale=1.2; 
% deltat = 0.0001 ; 
deltat = 0.08 ;
% deltat=.01 ;
geps=.001*h0; deps=sqrt(eps)*h0;
densityctrlfreq=30;

if isempty(p)
    %% 1. Create initial distribution in bounding box (equilateral triangles)
    [x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
    x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
    p=[x(:),y(:)];                                       % List of node coordinates

    %% 2. Remove points outside the region, apply the rejection method
    p=p(feval(fd,p,0)<geps,:);                      % Keep only d<0 points
    r0=1./feval(fh,p).^2;                           % Probability to keep point
    max_r0 = 1/h0^2;    %max_r0 = max(r0);
    p=p(rand(size(p,1),1) < r0/max_r0,:);           % Rejection method
end
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
pfix=unique(pfix,'rows'); nfix=size(pfix,1);
% return ; 
p=[pfix; p];                                         % Prepend fix points
N=size(p,1);                                         % Number of points N

%% Iterate
count=0; prcq_o = 0;
pold=inf;                                            % For first iteration
if plot_on == 1
    clf,view(2),axis equal,axis off
end
toc
while 1
  tic
  if mod(it,10) == 0
    disp(sprintf('Iteration = %d',it)) ;    
  end
  count=count+1;
  
  % 3. Retriangulation by the Delaunay algorithm
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    t=delaunayn(p);                                  % List of triangles
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
    t=t(feval(fd,pmid,0)<-geps,:);                   % Keep interior triangles
    
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    
    % 5. Graphical output of the current mesh
    if plot_on == 1
        cla,patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);
        drawnow
        if ~isempty(pfix)
            hold on
            movable = setdiff(p,pfix,'rows');
            plot(movable(:,1),movable(:,2),'o');
            hold off
        end
    end
    %
  end
  
  % Getting element quality and check goodness
  tq = gettrimeshquan( p, t);
  prcq = length(find(tq.qm > 0.6))/length(t);
  if it > itmax/4 && (prcq > 0.999 || (mod(it,50) == 0 && prcq <= prcq_o))
%       if any(tq.vang > 130/180*pi)
%           p(t(tq.vang > 130/180*pi),:) = [];
%           N=size(p,1); pold=inf;
%           continue;
%       end  
      endflag = remove_small_connectivity;
      if endflag == 0; continue; end 
      disp('quality of mesh is good enough, exit')
      break;
  end
  if mod(it,50) == 0
    prcq_o = prcq;
  end
  
  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
  hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2);
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  
  % Density control - remove points that are too close
  if mod(count,densityctrlfreq)==0   
      if any(L0>2*L)
        p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
        N=size(p,1); pold=inf;
        continue;
      end
  end
  
  F=max(L0-L,0);                                     % Bar forces (scalars)
  Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
  p=p+deltat*Ftot;                                   % Update node positions

  % 7. Bring outside points back to the boundary
  d=feval(fd,p,0); ix=d>0;                 % Find points outside (d>0)
  dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],0)-d(ix))/deps; % Numerical
  dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],0)-d(ix))/deps; % gradient 
  dgrad2=dgradx.^2+dgrady.^2;
  p(ix,:)=p(ix,:)-[d(ix).*dgradx./dgrad2,d(ix).*dgrady./dgrad2];    % Project

  % 8. Termination criterion: All interior nodes move less than dptol (scaled)
  it = it + 1 ;
  if ( it > itmax )
      endflag = remove_small_connectivity;
      if endflag == 0; continue; end 
      disp('too many iterations, exit')
      break ; 
  end
   
  if ( max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0) < dptol )
      endflag = remove_small_connectivity;
      if endflag == 0; continue; end 
      disp('no movement of nodes, exit')
      break; 
  end
  
  toc
end

%% Clean up and plot final mesh
[p,t] = fixmesh(p,t);

% DW, get rid of element with vanish area
xr = 0.5*(p(t(:,2),:) - p(t(:,1),:)) ; 
xs = 0.5*(p(t(:,3),:) - p(t(:,1),:)) ;

JJ = (xr(:,1).*xs(:,2) - xs(:,1).*xr(:,2)) ;
iflatel = find( JJ < 20*eps ) ;

inonflatel = setdiff(1:length(JJ),iflatel) ; 
t = t(inonflatel,:) ;

% WP Fix bad interior angle elements
[ vtov, nnv, ~, ~, ~ ] = NodeConnect2D( t );
tq = gettrimeshquan( p, t);
[I,~] = find(tq.vang < 30*pi/180 | tq.vang > 130*pi/180);
%       
for i = I'
    % Get the node belonging to largest angle
    [ang,j] = max(tq.vang(i,:));
    if ang < pi*90/180
        % Get node belonging to smallest angle
        [~,j] = min(tq.vang(i,:));
    end
    node = t(i,j) ;
    % Move this node to center of surrounding nodes
    center = mean(p(vtov(1:nnv(node),node),:));
    p(node,:) = center;
end
%
if plot_on == 1
    simpplot(p,t)
end

function endflag = remove_small_connectivity
    % Get node connectivity (look for 4)
    [ ~, ~, ~, nn, ~ ] = NodeConnect2D( t );
    % Make sure they are not boundary nodes
    d=feval(fd,p,0);
    if isempty(find(nn' == 4 & d < -h0,1))
        endflag = 1;
    else
        endflag = 0;
        p(nn' == 4 & d < -h0,:) = [];
        N=size(p,1); pold=inf;
    end
    return;
end

end