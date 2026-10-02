
clear all
addpath('./data')
addpath('./customcolormap/')
load('oscfoil_data.mat')
load('foilcoords.mat')
load("opdmd_modes.mat")
mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});

%% Dyn Sys Full Size Original
close all
width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = 0.5;
umax = 1.5;
xt = 2; %tip of the foil %1.4;
yt = max(y)/2;
SetFont = 'Times New Roman';
ps = 99-12; %start time where foil is flat
pnt = 99;
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
ind = [2, 4]
r=1;
for i=period_sq([2,4])%:800
    fig=figure
    set(gcf,'Position',[100 100 width height])
    uplot = u(:,:,i); %u2(:,:,i); %u1(:,:,i) + umean; %;
    contourf(xMat - xt,yMat - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(mycolormap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k';
    c.TickLength = .02;
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'w','EdgeColor','k')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'k'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k'; 
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['uorig_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i
end

%% Dyn Sys Full Size modes
close all
width =600; 
height = 250;
axislength = 20;
ep=0.01;
%umin = 0.5;
%umax = 1.5;
%umin = 0;
%umax = 2;
umin = -0.2;
umax = 0.2;
xt = 2; %1.4;
yt = max(y)/2;
SetFont = 'Times New Roman';
ps = 99-12;
pnt = 99;
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
ind = [2, 4]
r=1;
for i=period_sq([2,4])%:800
    fig=figure
    set(gcf,'Position',[100 100 width height])
    uplot = u2(:,:,i); %u2(:,:,i); %u1(:,:,i)+umean; %u1(:,:,i) + umean; %;
    contourf(xMat - xt,yMat - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(mycolormap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'w','EdgeColor','k')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'k'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k'; 
    %xlim([-0.7,x(end) - xt])
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['mode_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i
end

%% FTLE baseline system
clear all
addpath('./data')
addpath('./slanCM')
addpath('./cse_ftledata')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load('CSEtest_m1mn_tseries.mat','xPos','yPos','cseIntegral','deltaInfty','sigma_ftle','dt','tLength')

FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
CSE = (deltaInfty.*cseIntegral).^2;
lhs = log(CSE)/abs(2*tLength*dt);
zeta = (1/abs(tLength*dt))*(log(deltaInfty)+log(cseIntegral)-log(sigma_ftle));
xl=xPos(:,:,1);
yl=yPos(:,:,1);
%% FTLE baseline system
width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = 0;
umax = 1;
xt = 2;
yt = max(y)/2;
SetFont = 'Times New Roman';
cmap = slanCM('gnuplot2');
cmap = customcolormap(linspace(0,1,5),{'#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %blue
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %bluegreen
ps = 99-12;
pnt = 99;
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
ind = [2, 4]
r=1;
for i=period_sq(ind)%:99 %50% 100%:100
    fig=figure
    set(gcf,'Position',[100 100 width height])
    uplot = FTLE(:,:,i);
    contourf(xl - xt,yl - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(cmap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'k','EdgeColor','w')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'w'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k';
    %xticklabels([])
    %yticklabels([])
    %xlim([-0.7,x(end) - xt])
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['FTLEp_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i
    %img = getframe(2);
    %writeVideo(vd,img);  
end
%close(vd)


%% Mode Sensitivity
clear all
addpath('./data')
addpath('./cse_ftledata')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load('CSEtest_m1mn_tseries.mat','xPos','yPos','cseIntegral','deltaInfty','sigma_ftle','dt','tLength')

FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
CSE = (deltaInfty.*cseIntegral).^2;
lhs = log(CSE)/abs(2*tLength*dt);
zeta = (1/abs(tLength*dt))*(log(deltaInfty)+log(cseIntegral)-log(sigma_ftle));
xl=xPos(:,:,1);
yl=yPos(:,:,1);
%% Mode Sensitivity
width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = -1;
umax = 0;
xt = 2;
yt = max(y)/2;
SetFont = 'Times New Roman';
cmap = customcolormap(linspace(0,1,5),{'#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %blue
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %bluegreen
ps = 99-12;
pnt = 99;
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
ind = [2, 4]
r=1;
for i=period_sq([2,4])%:99 %50% 100%:100
    fig=figure
    set(gcf,'Position',[100 100 width height])
    uplot = lhs(:,:,i);
    contourf(xl - xt,yl - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(cmap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'k','EdgeColor','w')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'w'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k';
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['MS_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i
end
%close(vd)

%% Zeta
clear all
addpath('./data')
addpath('./slanCM/')
addpath('./cse_ftledata')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load('CSEtest_m1mn_tseries.mat','xPos','yPos','cseIntegral','deltaInfty','sigma_ftle','dt','tLength')

FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
CSE = (deltaInfty.*cseIntegral).^2;
lhs = log(CSE)/abs(2*tLength*dt);
zeta = (1/abs(tLength*dt))*(log(deltaInfty)+log(cseIntegral)-log(sigma_ftle));
xl=xPos(:,:,1);
yl=yPos(:,:,1);
%% Zeta
close all
width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = -1.2;
umax = -0.8;
xt = 2;
yt = max(y)/2;
SetFont = 'Times New Roman';
cmap = slanCM('gnuplot2');
cmap = customcolormap(linspace(0,1,5),{'#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %blue
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %blue
%cmap = customcolormap(linspace(0,1,6),{'#184ADC','#1FB8FF','#1fe9ff','#54ff82','#ffffff','#ffffff'}); %blue
%cmap = slanCM('plasma');
%vd = VideoWriter(['zeta_bluegreen_m1mn.avi']);
%vd.FrameRate = 15;
%vd.Quality = 95;
%open(vd)
ps = 99-12;
pnt = 99;
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
r=1
ind = [2, 4]
for i=period_sq([2,4])%:99 %50% 100%:100
    fig=figure
    set(gcf,'Position',[100 400 width height])
    uplot = zeta(:,:,i);
    contourf(xl - xt,yl - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(cmap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'k','EdgeColor','w')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'w'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k';
    %xticklabels([])
    %yticklabels([])
    %xlim([-0.7,x(end) - xt])
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['zeta_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i
    %img = getframe(2);
    %writeVideo(vd,img);  
end
%close(vd)

%% FTLE full system
clear all
addpath('./data')
addpath('./slanCM')
addpath('./cse_ftledata')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load('bFTLE_m1m2mn_tseries.mat')
%% FTLE full system
FTLE = (1/abs(tLength*dt))*log(sigma);
xl=xPos(:,:,1);
yl=yPos(:,:,1);

width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = 0;
umax = 1;
xt = 2;
yt = max(y)/2;
SetFont = 'Times New Roman';
cmap = slanCM('gnuplot2');
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %bluegreen
ps = 99-12;
pnt = 99;
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
r=1
ind = [2, 4]
for i=period_sq([2,4])%:99 %50% 100%:100
    fig=figure
    set(gcf,'Position',[100 100 width height])
    uplot = FTLE(:,:,i);
    contourf(xl - xt,yl - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(cmap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'k','EdgeColor','w')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'w'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k';
    %xticklabels([])
    %yticklabels([])
    %xlim([-0.7,x(end) - xt])
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['FTLEf_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i 
end

%% FTLE orig system
clear all
addpath('./data')
addpath('./slanCM')
addpath('./cse_ftledata')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load('bFTLEorig_tseries.mat')
%% FTLE full system
FTLE = (1/abs(tLength*dt))*log(sigma);
xl=xPos(:,:,1);
yl=yPos(:,:,1);

width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = 0;
umax = 1;
xt = 2;
yt = max(y)/2;
SetFont = 'Times New Roman';
cmap = slanCM('gnuplot2');
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %bluegreen
ps = 99-12;
pnt = 99;
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
r=1
ind = [2, 4]
for i=period_sq([2,4])%:99 %50% 100%:100
    fig=figure
    set(gcf,'Position',[100 100 width height])
    uplot = FTLE(:,:,i);
    contourf(xl - xt,yl - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(cmap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'k','EdgeColor','w')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'w'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k';
    %xticklabels([])
    %yticklabels([])
    %xlim([-0.7,x(end) - xt])
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    exportgraphics(fig, ['orig_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i 
end

%%
clear all
addpath('./data')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load("opdmdmodes.mat")
mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
%% Velocity Original
close all
width =600; 
height = 250;
axislength = 20;
ep=0.01;
umin = 1.5;
umax = 0.5;
xt = 2; %tip of the foil %1.4;
yt = max(y)/2;
SetFont = 'Times New Roman';
ps = 99-12; %start time where foil is flat
pnt = 99;
figure
set(gcf,'Position',[100 100 width height])
period = [ps, ps + round(pnt/4), ps + round(pnt/2), ps + round(3*pnt/4) ]
period_sq = [ps + round(3*pnt/4), ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) ]
period_sq = [ps + round(3*pnt/4) + round(pnt/8), ps + round(3*pnt/4) + round(pnt/4), ps + round(3*pnt/4) + round(3*pnt/8), ps + round(3*pnt/4) + round(pnt/2) - pnt]-pnt
for i=1:400%:800
    w = curl(x,y,u(:,:,i),v(:,:,i));
    uplot = u(:,:,i); %u2(:,:,i); %u1(:,:,i) + umean; %;
    contourf(xMat - xt,yMat - yt,uplot,80,'LineStyle','none')
    colorbar
    colormap(mycolormap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.Color = 'k'
    c.TickLength = .02;
    hold on
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    fill([xfoil(foilStartpx:foilEndpx) ,foilEnd,foilStart] - xt, [yfoil(foilStartpx:foilEndpx,i+phase)',0,0]-yt,'w','EdgeColor','k')
    set(gcf,'color','w')
    ax2 = gca; ax2.FontSize = axislength; ax2.Color = 'k'; ax2.FontName = SetFont;
    ax2.XColor = 'k'; ax2.YColor = 'k'; 
    xlim([-1.3,x(end) - xt])
    hold off
    drawnow
    %exportgraphics(fig, ['upart_',num2str(ind(r)),'.tif'], 'Resolution', 1200)
    r = r+1;
    i
end

