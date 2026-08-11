%% MS Field with Phases Frames
%C:\Users\Gator\Documents\PhD Research 2023\JFM 2024\scripts\current
clear all
addpath('data\')
q=load("currentMSEsnapshot_4phaseonlyfix.mat");
p1 = load('currentboundmse_tff_location1_b4_T24.mat');
p2 = load('currentboundmse_tff_location2_b4_T24.mat');
 %%
 close all
lhs = log(q.MSE)/(2*q.T); %MSE Field
ep=0.01
omega = q.omega;
f = omega/(2*pi);
cmin = -0.3;
cmax = 0.2;
%cmin = 0;
%cmax = 0.25;
%cmin = -0.10;
%cmax = 0.15;
dt = q.dt;
Tosc = 1/f;
frameperiod = Tosc/dt;
axislength = 30
SetFont = 'Times New Roman'

close all
addpath('customcolormap\')
%load('ForwardFTLEgyreFlow.mat')
if q.int == 'f'
    %mycolormap = customcolormap(linspace(0,1,8),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff'});
    mycolormap = customcolormap(linspace(0,1,10),{'#FADA0E','#EE2020','#4A4A4A','#6E6E6E','#9A9A9A','#ffffff','#ffffff','#ffffff','#ffffff','#ffffff'});
    %mycolormap = customcolormap(linspace(0,1,8),{'#f5d586','#edb121','#ca6827','#b8412a','#a9212d','#000000','#000000','#000000'});
    %mycolormap = customcolormap(linspace(0,1,7),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff'});
else
    mycolormap = customcolormap(linspace(0,1,6),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff'});
end
mycolormap = customcolormap(linspace(0,1,5), flip({'#0d0780','#6e13a2','#bc537b','#ea9b51','#f1f758'})); %plasma (for zeta pub)
mycolormap = customcolormap(linspace(0,1,7),{'#FADA0E','#EE2020','#000000','#000000','#000000','#000000','#000000'}); %red

px1 = p1.xgmin;
py1 = p1.ygmin;
px2 = p2.xgmin;
py2 = p2.ygmin;
for i=[2,51,101,151]%[1, frameperiod/4, frameperiod/2, 3*frameperiod/4]%:size(q.MSE,3)
    fig=figure
    set(gcf,'Position',[100 100 1000 500])
    i
    contourf(q.x0grid',q.y0grid',lhs(:,:,i),100,'LineStyle','none')
    hold on
    scatter(px1, py1, 100,'filled','g'); % 'filled' just makes the markers solid
    scatter(px2, py2, 100,'filled','g'); % 'filled' just makes the markers solid
    axis([q.xgmin q.xgmax q.ygmin q.ygmax])
    %yticks([])
    %xticks([])
    axis equal
    yticks([0, 0.5, 1.0])
    xticks([0, 0.5, 1.0, 1.5, 2.0])
    yticklabels({'0', '0.5', '1.0'})
    xticklabels({'0', '0.5', '1.0', '1.5', '2.0'})
    axis([q.xgmin q.xgmax q.ygmin q.ygmax])
    ax2 = gca; ax2.FontName = SetFont; ax2.FontSize = axislength
    %colormap([flip(spring);flip(cool)])
    %colormap(jet)
    set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
    c=colorbar
    c.FontSize =axislength;
    c.FontName = 'Times New Roman'
    colormap(mycolormap)
    clim([cmin,cmax]) %zeta
    c.Ticks = linspace(cmin,cmax,3);
    %c.Ticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0];
    set(gcf,'color','w')
    c.LineWidth = 1.2;
    c.TickLength = .02;
    set(gcf,'color','w')
    %axis off
    drawnow
    %exportgraphics(fig, ['msfield_location1_t',num2str(i),'.tif'], 'Resolution', 1200)
end


