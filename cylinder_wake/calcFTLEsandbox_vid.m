%% DMD Recon Flow Field
% C:\Users\Gator\Documents\PhD Research 2023\Lagrangian Coherent Structures\Haller_2020\experimental data
clear all
addpath('./data')
addpath('./customcolormap/')
load('cylindervelocity_zen.mat')
%load("cylinderPODmodesPair.mat")
%%
ufq = u;
vfq = v;
clear ufint vfint ugint vgint u1 u2 u3 u4 umean v1 v2 v3 v3 v4 vmean xnMat ynMat
%%
uMesh=ufq;
vMesh=vfq;
clear ufq vfq ugq vgq

%% Compute FTLE
clear sigmaS

xMinROI = -2;%-5;
xMaxROI = 20;
yMinROI = -5;%-7;
yMaxROI = 5;%7;
ROInx = 800;
ROIny = 400;
method = 'RK4';%'RK4';
ii=1;
tVec = (0:1:(nt-1))*dt';
yVec=y; %y(1,:)';
xVec=x;%x(1,256:769)';
%%
% Extrapolation funtions for FTLE
Fx = interpolantxy(uMesh,xVec,yVec,tVec,1);
Fy = interpolantxy(vMesh,xVec,yVec,tVec,0);

uExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx);
vExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy);
%%
clc
%Forwards FTLE
%inteval of integration from t0 to t
 % tLength = 240; %two periods %negative if backwards int
 % finc = 1; 
 % tStep = 1; %negative if backwards int
 % fstart = 120+1; %240+1;
 % %Last Frame to compute FTLE Field
 % fend = 120+120+1; %240+120+1;

%Backwards FTLE
tLength = -199; %negative if backwards int
finc = 1; %negative if backwards int
tStep = -1;
fstart = 200; %105;
%Last Frame to compute FTLE Field
fend = "none"; %105+60; %4*pnt;
fLoop = "none"; %fstart:finc:fend;
sigma = zeros([ROInx, ROIny]);
save('bFTLE_original.mat','sigma','tLength','tVec','dt','xMaxROI','xMinROI','yMaxROI','yMinROI','fstart','method',"-append")
%%%%%%%%%%%
%%
for tStart = fstart %fLoop
    [sigma, xPos, yPos] = FTLEsb(uMesh, vMesh, ...
    xVec, yVec, ...
    tStart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ...
    ROInx, ROIny, method, ...
    'extrap', true, 'uExtrap', uExtrap, 'vExtrap', vExtrap);
    %save('fFTLE_ss_orig.mat','sigma','xPos','yPos','-append')
end

save('bFTLE_original.mat','sigma','xPos','yPos')

%save('bFTLE_ss_orig.mat','sigma','xPos','yPos','tLength','tVec','dtn','fstart','finc','fend')
%% FTLE only
Tperiod = (abs(tLength)*dtn);
FTLE = (1/Tperiod)*log(sigma);
xg=xPos(:,:,1);
yg=yPos(:,:,1);
% FTLE
mycolormap = customcolormap(linspace(0,1,7),{'#FADA0E','#EE2020','#000000','#000000','#000000','#000000','#000000'}); %red
mycolormap = customcolormap(linspace(0,1,7),{'#62dcfc','#184ADC','#10349e','#0c297d','#081f5e','#000000','#000000'});
width =800; 
height = 400;
axislength = 20;
ep=0.01;
umin = 0;
umax = 0.3;
%umin = 0.01;
%umax = 0.04;
SetFont = 'Times New Roman';
figure(1)
set(gcf,'position',[200,100,width,height])
%vd = VideoWriter(['fFTLE_original.avi']);
%vd.FrameRate = 12;
%open(vd)
for i=1:120 %50% 100%:100
    FTLEa = FTLE(:,:,i);
    %FTLEa(FTLEa<0.2)=0;
    contourf(xg,yg,FTLEa,80,'LineStyle','none')
    xlim([-1,8])
    %xlim([-5,4])
    colorbar
    colormap(mycolormap)
    %colormap(jet)
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
    fill(xcyl,ycyl,[0.4 0.4 0.4])  % place cylinder
    %plot(xcyl,ycyl,'w','LineWidth',1.4) % cylinder boundary
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)
%print('-vector','-depsc','template.eps')
%print('-vector','-dtiffn','ftle_backwards_original.tif')
%%
vd = VideoWriter(['test.avi']);
vd.FrameRate = 15;
open(vd)
figure(1)
set(1,'position',[10 100 1400 1000])
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
frames=fstart:1:fstart+tLength
r=1;
for i=frames
    i
    %vort = curl(xqMat(:,:,1),yqMat(:,:,1),uq(:,:,i),vq(:,:,i));
    contourf(xqMat(:,:,1),yqMat(:,:,1),uq(:,:,i),200,'linestyle','none')
    colormap(mycolormap)
    c=colorbar;
    %clim([-1,1])
    clim([0,1.2])
    hold on
    step=4;
    %ind = frames(1)-frames(fstart-i+1)+1;
    scatter(xPos(1:step:end,1:step:end,r),yPos(1:step:end,1:step:end,r),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerFaceAlpha',0.1)
    axis equal
    set(gcf,'color','white')
    %xlim([xMinROI xMaxROI + 5])
    %ylim([yMinROI-0.5 yMaxROI+0.5])
    xlim([xMinROI xMaxROI+5])
    ylim([yMinROI-0.2 yMaxROI+0.2])
    hold off
    drawnow
    r=r+1;
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

