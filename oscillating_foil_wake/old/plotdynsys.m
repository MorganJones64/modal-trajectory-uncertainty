
clear all
addpath('./data')
addpath('./fxn')
addpath('./customcolormap/')
load('oscfoilfulldata.mat')
load('foilcoords.mat')
load("opdmdmodes.mat")
mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
%%
width =800; 
height = 300;
axislength = 20;
ep=0.01;
%umin = 0;
%umax = 2;
umin = -1;
umax = 1;
SetFont = 'Times New Roman';
figure(1)
set(gcf,'Position',[100 100 width height])
%vd = VideoWriter(['velocityvert_original.avi']);
%vd.FrameRate = 15;
%open(vd)
for i=1%:151 %50% 100%:100
    uplot = vmean(:,:,i)+v1(:,:,i)+v2(:,:,i)+v3(:,:,i)+v4(:,:,i);
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
    fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,i+phase)',0,0],'w')
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)

%%
width =800; 
height = 300;
axislength = 20;
ep=0.01;
umin = 0;
umax = 2;
%umin = -1;
%umax = 1;
SetFont = 'Times New Roman';

figure(2)
set(gcf,'Position',[100 100 width height])
%vd = VideoWriter(['velocityvert_original.avi']);
%vd.FrameRate = 15;
%open(vd)
for i=1:151 %50% 100%:100
    uplot = u(:,:,i);
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
    fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,i+phase)',0,0],'w')
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    xlim([0.8,x(end)])
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)


%%
width =800; 
height = 400;
axislength = 20;
ep=0.01;
umin = 0;
umax = 2;
umin = -2;
umax = 2;
SetFont = 'Times New Roman';

figure(2)
set(gcf,'Position',[100 100 width height])
%vd = VideoWriter(['velocityvert_original.avi']);
%vd.FrameRate = 15;
%open(vd)
for i=1:151 %50% 100%:100
    uplot = u(:,:,i);
    vort = curl(xMat(:,:,1),yMat(:,:,1),u(:,:,i),(-1)*v(:,:,i));
    contourf(xMat(:,:,1),yMat(:,:,1),vort,200,'linestyle','none')
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
    fill([xfoil(foilStartpx:foilEndpx),foilEnd,foilStart], [yfoil(foilStartpx:foilEndpx,i+phase)',0,0],'w')
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)






%%
clear all
addpath('./data')
addpath('./customcolormap/')
load('cylindervelocity.mat')
load('cylinderDMDmodes.mat')
mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
%%
width =800; 
height = 400;
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
    uplot = umean+u1(:,:,i)+u2(:,:,i)+u3(:,:,i)+u4;
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
    fill(xcyl,ycyl,[.6 .6 .6])  % place cylinder
    plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
    ax2 = gca; axis equal; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
    hold off
    drawnow
    i
    %img = getframe(1);
    %writeVideo(vd,img);  
end
%close(vd)