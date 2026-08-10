
clear all
addpath('./data')
addpath('./fxn')
addpath('./customcolormap/')
load('cylindervelocity_zen.mat')
load('cylinderPODmodesPair.mat')
load('cylcoords.mat')
%%
mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
%mycolormap = customcolormap(linspace(0,1,12),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
%mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});

width =800; 
height = 320;
axislength = 20;
ep=0.01;
umin = -0.04;
umax = 0.04;
% umin = 0;
% umax = 1.2;
SetFont = 'Times New Roman';
figure(1)
set(gcf,'Position',[100 100 width height])
%vd = VideoWriter(['velocityvert_original.avi']);
%vd.FrameRate = 15;
%open(vd)
for i=1%:151 %50% 100%:100
    uplot = u3(:,:,end);
    contourf(xMat,yMat,uplot,100,'LineStyle','none')
    colorbar
    colormap(mycolormap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    xlabel('$x$','Interpreter','latex','fontsize',axislength);
    ylabel('$y$','Interpreter','latex','fontsize',axislength);
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    set(gcf,'color','w')
    fill(xcyl,ycyl,[.3 .3 .3])  % place cylinder
    plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
    xlim([-2,20])
    ylim([-4,4])
    xticks(0:4:20)
    yticks(-4:4:4)
    ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)

%%
width =800; 
height = 320;
axislength = 20;
ep=0.01;
umin = 0;
umax = 1.2;
SetFont = 'Times New Roman';
figure(1)
set(gcf,'Position',[100 100 width height])
%vd = VideoWriter(['velocityvert_original.avi']);
%vd.FrameRate = 15;
%open(vd)
for i=1%:151 %50% 100%:100
    %uplot = u(:,:,i);%umean+u1(:,:,i)+u2(:,:,i)+u3(:,:,i)+u4(:,:,i);
    %contourf(xMat,yMat,uplot,80,'LineStyle','none')
    uplot = umean+u1(:,:,end)+u2(:,:,end);%umean+u1(:,:,i)+u2(:,:,i)+u3(:,:,i)+u4(:,:,i);
    contourf(xMat,yMat,uplot,80,'LineStyle','none')
    colorbar
    colormap(mycolormap)
    axis equal
    clim([umin,umax])
    c=colorbar;
    c.Ticks =linspace(umin,umax,3);
    c.FontSize =axislength;
    c.FontName = 'Times New Roman';
    c.LineWidth = 1.2;
    c.TickLength = .02;
    %xticks([3,4,5,6,7,8,9,10,11,12])
    hold on
    xlabel('$x$','Interpreter','latex','fontsize',axislength);
    ylabel('$y$','Interpreter','latex','fontsize',axislength);
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    set(gcf,'color','w')
    xlim([-2,20])
    ylim([-4,4])
    xticks(0:4:20)
    yticks(-4:4:4)
    fill(xcyl,ycyl,[.6 .6 .6])  % place cylinder
    plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
    ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)