%cd(fileparts(matlab.desktop.editor.getActiveFilename))
%% Get DAta
clear all
addpath('./data')
addpath('./fxn')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load("opdmdmodes.mat")
%% Compute FTLE
clear sigmaS

xMinROI = min(x);
xMaxROI = max(x);
yMinROI = min(y);
yMaxROI = max(y);
ROInx = 200;
ROIny = 100;
method = 'RK4';%'RK4';
ii=1;

%uMesh = umean+u1+u2+u3+u4;
%vMesh = vmean+v1+v2+v3+v4;
uMesh = u;
vMesh = v;
umeanprof = linspace(1,1,ny); %mean(uMesh(:,end,:),3)'; %
vmeanprof = linspace(0,0,ny); %mean(vMesh(:,end,:),3)'; %
yVec=y';
xVec=x';
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
% tLength = 60; %300; %round(6*pnt); 
% tStep = 1;
% finc = 1;
% %First Frame to compute FTLE Field
% fstart = 1;
% %Last Frame to compute FTLE Field
% fend = 1; %4*pnt;
pnt = 99;

%Backwards FTLE
 tLength = -round(0.5*pnt); %negative if backwards int
 finc = 1; 
 tStep = -1; %negative if backwards int
 fstart = 2*pnt;
 %Last Frame to compute FTLE Field
 fend = 3*pnt+1; %4*pnt;

% Frame Increment
% Frame loop vector
fLoop = fstart:finc:fend;
% Preallocate space for sigma
sigma = zeros([ROInx, ROIny, length(fLoop)]);
%%%%%%%%%%%
%%
for tStart = fLoop
((tStart-fstart)/(fend-fstart))*100
[sigma(:, :, tStart-fLoop(1)+1), xPos, yPos] = FTLEsb(uMesh, vMesh, ...
    xVec, yVec, ...
    tStart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ...
    ROInx, ROIny, method, ...
    'extrap', true, 'uExtrap', uExtrap, 'vExtrap', vExtrap);
end

%save('bFTLE_ss.mat','sigma','xPos','yPos','tLength','tVec','dt','fstart','finc','fend')
%% FTLE only
Tperiod = (abs(tLength)*dt);
FTLE = (1/Tperiod)*log(sigma);
xg=xPos(:,:,1);
yg=yPos(:,:,1);
%%
% FTLE
load('foilcoords.mat')
mycolormap = customcolormap(linspace(0,1,10),{'#FADA0E','#EE2020','#4A4A4A','#6E6E6E','#9A9A9A','#ffffff','#ffffff','#ffffff','#ffffff','#ffffff'});
mycolormap = customcolormap(linspace(0,1,11),{'#1FB8FF','#184ADC','#10349e','#4A4A4A','#6E6E6E','#9A9A9A','#BDBDBD','#ffffff','#ffffff','#ffffff','#ffffff'});
mycolormap = customcolormap(linspace(0,1,11),{'#8ad1e1','#184ADC','#10349e','#4A4A4A','#6E6E6E','#9A9A9A','#BDBDBD','#ffffff','#ffffff','#ffffff','#ffffff'});
mycolormap = customcolormap(linspace(0,1,11),{'#56c1de','#184ADC','#10349e','#4A4A4A','#6E6E6E','#9A9A9A','#BDBDBD','#ffffff','#ffffff','#ffffff','#ffffff'});
mycolormap = customcolormap(linspace(0,1,8),{'#56c1de','#184ADC','#10349e','#4A4A4A','#6E6E6E','#9A9A9A','#BDBDBD','#ffffff'});
%mycolormap = customcolormap(linspace(0,1,6),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff'})
%mycolormap = customcolormap(linspace(0,1,7),flip({'#000000','#000000','#a9212d','#b8412a','#ca6827','#edb121','#f5d586'}));
%mycolormap = customcolormap(linspace(0,1,5),{'#FADA0E','#EE2020','#4A4A4A','#9A9A9A','#ffffff'});
%mycolormap = customcolormap(linspace(0,1,8),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff','#ffffff','#ffffff'});
%mycolormap = customcolormap(linspace(0,1,8),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff'});
%mycolormap = customcolormap(linspace(0,1,8),{'#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d','#000000','#000000','#000000'});
%mycolormap = customcolormap(linspace(0,1,8),{'#a9d9f8','#2295e0','#0758ab','#0d449d','#000000','#000000','#000000','#000000'});
%mycolormap = customcolormap(linspace(0,1,8),{'#a9d9f8','#2295e0','#0758ab','#000000','#000000','#000000','#000000','#000000'});
width =800; 
height = 400;
axislength = 20;
ep=0.01;
umin = 0;
umax = 1;
%umin = 0.01;
%umax = 0.04;
SetFont = 'Times New Roman';
figure
set(gcf,'position',[200,100,width,height])
%v = VideoWriter(['bFTLE_original.avi']);
%v.FrameRate = 10;
%open(v)
for i=1%:100 %50% 100%:100
    FTLEa = FTLE(:,:,i);
    %FTLEa(FTLEa<0.2)=0;
    contourf(xg,yg,FTLEa,80,'LineStyle','none')
    axis([xMinROI xMaxROI yMinROI yMaxROI])
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
    fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,fstart+phase)',0,0],'w')
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    drawnow
    i
    %img = getframe(1);
    %writeVideo(v,img);  
end
%close(v)
%print('-vector','-depsc','template.eps')
%print('-vector','-dtiffn','ftle_backwards_original.tif')



%%
% FTLE
load('foilcoords.mat')
mycolormap = customcolormap(linspace(0,1,8),{'#56c1de','#184ADC','#10349e','#4A4A4A','#6E6E6E','#9A9A9A','#BDBDBD','#ffffff'});

width =800; 
height = 400;
axislength = 20;
ep=0.01;
umin = 0;
umax = 1;
%umin = 0.01;
%umax = 0.04;
SetFont = 'Times New Roman';
figure(1)
set(gcf,'position',[200,100,width,height])
v = VideoWriter(['bFTLE_original.avi']);
v.FrameRate = 10;
open(v)
for i=1:99 %50% 100%:100
    FTLEa = FTLE(:,:,i);
    %FTLEa(FTLEa<0.2)=0;
    contourf(xg,yg,FTLEa,80,'LineStyle','none')
    axis([xMinROI xMaxROI yMinROI yMaxROI])
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
    fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,fLoop(i)+phase)',0,0],'w')
    hold off
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    drawnow
    i
    img = getframe(1);
    writeVideo(v,img);  
end
close(v)
%print('-vector','-depsc','template.eps')
%print('-vector','-dtiffn','ftle_backwards_original.tif')
%%
vd = VideoWriter(['bparticleAdvection2.avi']);
vd.FrameRate = 20;
open(vd)
figure(1)
set(1,'position',[10 100 1400 1000])
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
frames=fstart:-1:fstart+tLength
for i=frames
    i
    vort = curl(xMat(:,:,1),yMat(:,:,1),u(:,:,i),v(:,:,i));
    contourf(xMat(:,:,1),yMat(:,:,1),vort,200,'linestyle','none')
    colormap(mycolormap)
    c=colorbar;
    clim([-1,1])
    step=4;
    ind = frames(1)-frames(fstart-i+1)+1;
    hold on
    scatter(xPos(1:step:end,1:step:end,ind),yPos(1:step:end,1:step:end,ind),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerFaceAlpha',0.2)
    axis equal
    set(gcf,'color','white')
    xlim([xMinROI-5 xMaxROI])
    ylim([yMinROI-0.2 yMaxROI+0.2])
    hold off
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

