
%% DMD Recon Flow Field
% C:\Users\Gator\Documents\PhD Research 2023\Lagrangian Coherent Structures\Haller_2020\experimental data
clear all
addpath('./data')
addpath('./fxn')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load("opdmdmodes.mat")
%%
ufMesh=umean+u1+u2;
vfMesh=vmean+v1+v2;
ugMesh=u3;
vgMesh=v3;
%[ny,nx,nt]=size(ur);
clear ufq vfq ugq vgq

%% Compute FTLE
clear sigmaS

xMinROI = min(x);
xMaxROI = max(x);
yMinROI = min(y);
yMaxROI = max(y);
ROInx = 300;
ROIny = 200;
method = 'RK4';%'RK4';
ii=1;

umeanprof = linspace(1,1,ny); %mean(uMesh(:,end,:),3)'; %
vmeanprof = linspace(0,0,ny); %mean(vMesh(:,end,:),3)'; %
yVec=y';
xVec=x';
%%
Ffx = interpolantxy(ufMesh,xVec,yVec,tVec,1);
Ffy = interpolantxy(vfMesh,xVec,yVec,tVec,0);
Fgx = interpolantxy(ugMesh,xVec,yVec,tVec,0);
Fgy = interpolantxy(vgMesh,xVec,yVec,tVec,0);

ufExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Ffx);
vfExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Ffy);
ugExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fgx);
vgExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fgy);
clear Ffx Ffy Fgx Fgy
%%
%%
clc
%Forwards FTLE
%inteval of integration from t0 to t
% tLength = 60; %300; %round(6*pnt); 
% tStep = 1;
% finc = 1;
% %First Frame to compute FTLE Field
% fstart = 1;
% %Last Frame to compute FTLE Field
% fend = 1; %4*pnt;
pnt = 99;
%Backwards FTLE
 tLength = -round(1.5*pnt); %negative if backwards int
 finc = 1; 
 tStep = -1; %negative if backwards int
 fstart = 2*pnt;
 %Last Frame to compute FTLE Field
 fend = 3*pnt+1; %4*pnt;
 fLoop = fstart:finc:fend;
%%
% Preallocate space for sigma
sigma_ftle = zeros([ROInx, ROIny, length(fLoop)]);
cseIntegral = zeros([ROInx, ROIny, length(fLoop)]);
deltaInfty = zeros([ROInx, ROIny, length(fLoop)]);
notes = 'f=m1+m2+mean, g=m3, backward ftle'
%%%%%%%%%%%
%save('CSEtest_m1m2mn.mat','sigma_ftle','cseIntegral','deltaInfty','dt','tLength','xMaxROI','xMinROI','yMaxROI','yMinROI','fstart','fend','fLoop','method','notes')
%% Continue
%fLoop = t0:finc:fend;
%%
for t0 = fstart %fLoop %t0 = tStart
    tic
    [sigma_ftle(:, :, t0-fLoop(1)+1), cseIntegral(:, :, t0-fLoop(1)+1), deltaInfty(:, :, t0-fLoop(1)+1), xPos, yPos] = CSEsb(ufMesh, vfMesh, ...
        ugMesh, vgMesh,...
        xVec, yVec, ...
        t0, tLength, tStep, fstart, eps, dt, ...
        xMinROI, xMaxROI, yMinROI, yMaxROI, ...
        ROInx, ROIny, method, ...
        'extrap', true, 'ufExtrap', ufExtrap, 'vfExtrap', vfExtrap, ...
        'ugExtrap', ugExtrap, 'vgExtrap', vgExtrap);
    toc
    save('CSEtest_m1m2mn.mat','sigma_ftle','cseIntegral','deltaInfty','xPos','yPos',"-append")
end

%%
clear all
addpath('csedata\')

%%
FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
CSE = (deltaInfty.*cseIntegral).^2;
lhs = log(CSE)/abs(2*tLength*dt);
zeta = (1/abs(tLength*dt))*(log(deltaInfty)+log(cseIntegral)-log(sigma_ftle));
zeta2 = log(deltaInfty.*cseIntegral./sigma_ftle);
rhs = FTLE + zeta;
x=xPos(:,:,1);
y=yPos(:,:,1);
%%
% Plot the MSE Field with FTLE field
load('foilcoords.mat')
addpath('customcolormap\')
%if q.int == 'f'
    %mycolormap = customcolormap(linspace(0,1,6),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff'});
%else
    mycolormap = customcolormap(linspace(0,1,6),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff'});
%end
%close all
ep=0.01;
SetFont = 'Times New Roman';
width =800; 
height = 400;
axislength = 20;

%umin = -0.5;
%umax = 0.2;

%umin = -0.6;
%umax = -0.35;

%umin = -0.6;
%umax = -0.5;

umin = 0;
umax = 1.5;

figure
set(gcf,'Position',[100 100 width height])
contourf(x,y,FTLE(:,:,1),200,'LineStyle','none')
axis([xMinROI xMaxROI yMinROI yMaxROI])
colorbar
%colormap(mycolormap)
colormap(jet)
axis equal
clim([umin,umax])
c=colorbar;
c.Ticks =linspace(umin,umax,3);
c.FontSize =axislength;
c.FontName = 'Times New Roman';
c.LineWidth = 1.2;
c.TickLength = .02;
hold on
fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,fstart+phase)',0,0],'w')
xlabel('$x$','Interpreter','latex','fontsize',axislength);
ylabel('$y$','Interpreter','latex','fontsize',axislength);
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
%fill(xcyl,ycyl,[1 1 1])  % place cylinder
%plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
drawnow

%% CSE
ep=0.01;
SetFont = 'Times New Roman';
width =800; 
height = 400;
axislength = 20;

%umin = -0.3;
%umax = -0.18;

umin = -0.3;
umax = -0.1;

umin = -1.2;
umax = -0.8;

%umin = 0;
%umax = 0.2;

figure
set(gcf,'Position',[100 100 width height])
contourf(x,y,zeta(:,:,1),200,'LineStyle','none')
axis([xMinROI xMaxROI yMinROI yMaxROI])
colorbar
%colormap(mycolormap)
colormap(jet)
axis equal
clim([umin,umax])
c=colorbar;
c.Ticks =linspace(umin,umax,3);
c.FontSize =axislength;
c.FontName = 'Times New Roman';
c.LineWidth = 1.2;
c.TickLength = .02;
hold on
xlabel('$x$','Interpreter','latex','fontsize',axislength);
ylabel('$y$','Interpreter','latex','fontsize',axislength);
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,fstart+phase)',0,0],'w')
ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
drawnow


%% CSE
ep=0.01;
SetFont = 'Times New Roman';
width =800; 
height = 400;
axislength = 20;

umin = -0.5
umax = 0.5;
figure
set(gcf,'Position',[100 100 width height])
contourf(x,y,lhs(:,:,1),200,'LineStyle','none')
axis([xMinROI xMaxROI yMinROI yMaxROI])
colorbar
%colormap(mycolormap)
colormap(jet)
axis equal
clim([umin,umax])
c=colorbar;
c.Ticks =linspace(umin,umax,3);
c.FontSize =axislength;
c.FontName = 'Times New Roman';
c.LineWidth = 1.2;
c.TickLength = .02;
hold on
xlabel('$x$','Interpreter','latex','fontsize',axislength);
ylabel('$y$','Interpreter','latex','fontsize',axislength);
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,fstart+phase)',0,0],'w')
ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
drawnow




%%
vd = VideoWriter(['bparticleAdvectionfps.avi']);
vd.FrameRate = 10;
open(vd)
figure(1)
set(1,'position',[10 100 1400 400])
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
frames=151:-1:1
for i=frames
    i
    %vort = curl(xMat,yMat,u(:,:,i),v(:,:,i));
    %contourf(xMat,yMat,vort,200,'linestyle','none')
    %colormap(mycolormap)
    %c=colorbar
    %clim([-1,1])
    %clim([0,1.5])
    %hold on
    step=4;
    ind = frames(1)-frames(151-i+1)+1
    scatter(xPos(1:step:end,1:step:end,ind),yPos(1:step:end,1:step:end,ind),'o','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerFaceAlpha',0.2)
    axis equal
    set(gcf,'color','white')
    %xlim([xMinROI xMaxROI + 5])
    %ylim([yMinROI-0.5 yMaxROI+0.5])
    xlim([xMinROI-5 xMaxROI])
    ylim([yMinROI-0.2 yMaxROI+0.2])
    %hold off
    drawnow
    img = getframe(1);
    writeVideo(vd,img);  
end
close(vd)

%%
function F = interpolantxy(uvec,xVec,yVec,tVec,bcs)
    nt = length(tVec);
    nx = length(xVec);
    ny = length(yVec);
    yb = 1;
    xb = 1;
    xVecExt = [xVec(1)-xb, xVec];
    yVecExt = [yVec(1)-yb, yVec, yVec(end)+yb];
    if bcs == 0
        uleft = zeros(ny,1);
        utb = zeros(1,nx+1);
    else
        uleft = ones(ny,1);
        utb = ones(1,nx+1);
    end
    for i = 1:nt
        utemp = [uleft, uvec(:,:,i)]; %fix
        uExt(:,:,i) = [utb; utemp; utb];
    end
    [yi,xi,ti] = ndgrid(yVecExt,xVecExt, tVec);
    F = griddedInterpolant(yi, xi, ti, uExt,'makima','nearest'); 
end

function u = uExtrapolateM(x, y, t, u, F) 
    u = F(y,x,t);
end

