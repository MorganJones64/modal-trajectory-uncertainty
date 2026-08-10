%%
clear all
addpath('./data')
addpath('./customcolormap')
load("CSE_fm1_2_gm3_DNS_800_400.mat")
load("cylcoords.mat")
%%
FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
CSE = (deltaInfty.*cseIntegral).^2;
lhs = log(CSE)/abs(2*tLength*dt);
zeta = (1/abs(tLength*dt))*(log(deltaInfty)+log(cseIntegral)-log(sigma_ftle));
zeta2 = log(deltaInfty.*cseIntegral./sigma_ftle);
rhs = FTLE + zeta;
x=xPos(:,:,1);
y=yPos(:,:,1);
%% CSE full
ep=0.01;
SetFont = 'Times New Roman';
width =800; 
height = 320;
axislength = 20;
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %bluegreen
cmap = customcolormap(linspace(0,1,5),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000'}); %bluegreen
field = zeta;
%field(FTLE<0.10) = nan;
% Example point coordinates
% px = [px1, px2, px3, px4]; % all x-coordinates
% py = [py1, py2, py3, py4]; % all y-coordinates
% 
% px = [2.02, 4.0, 6.00, 6.37, px9, px10, px11];
% py = [-0.41, -0.44, -0.589, 0.138, py9, py10, py11];

floc = './data/boundextra/';
floc2 = './data/';
addpath('./data')
addpath('./data/boundextra/')

p61 = load([floc,'particlePos_backwards_p13_original.mat']);

px6 = p61.xMinROI + p61.dtol
py6 = p61.yMinROI + p61.dtol 

p51 = load([floc,'particlePos_backwards_p12_original.mat']);

px5 = p51.xMinROI + p51.dtol
py5 = p51.yMinROI + p51.dtol 

p41 = load([floc,'particlePos_backwards_p11_original.mat']);

px4 = p41.xMinROI + p41.dtol
py4 = p41.yMinROI + p41.dtol 

p31 = load([floc,'particlePos_backwards_p10_original.mat']);

px3 = p31.xMinROI + p31.dtol
py3 = p31.yMinROI + p31.dtol

p21 = load([floc2,'particlePos_backwards_p6c_original.mat']);

px2 = p21.xMinROI + p21.dtol
py2 = p21.yMinROI + p21.dtol

p11 = load([floc2,'particlePos_backwards_p5_original.mat']);

px1 = p11.xMinROI + p11.dtol
py1 = p11.yMinROI + p11.dtol 

px = [px1 px2 px3 px4 px5 px6];
py = [py1 py2 py3 py4 py5 py6];

umin = -0.45;
umax = -0.25;

umin = -0.4;
umax = -0.25;

figure
set(gcf,'Position',[100 100 width height])
contourf(x,y,field(:,:,1),200,'LineStyle','none')
hold on
scatter(px, py, 20,'filled','MarkerFaceColor',[255, 0, 0]/255,'MarkerFaceAlpha',0.7); % 'filled' just makes the markers solid
axis([xMinROI xMaxROI yMinROI yMaxROI])
colorbar
%colormap(mycolormap)
colormap(cmap)
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
fill(xcyl,ycyl,[1 1 1])  % place cylinder
plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
xlim([-2,20])
ylim([-4,4])
xticks(0:4:20)
yticks(-4:4:4)
%xlim([-1,8])
%ylim([-2,2])
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
drawnow

%% CSE w cropped FTLE field
ep=0.01;
SetFont = 'Times New Roman';
width =800; 
height = 320;
axislength = 20;
cmap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'}); %bluegreen
field = zeta;
field(FTLE<0.12) = nan;
% Example point coordinates
px = [px1, px2, px3, px4]; % all x-coordinates
py = [py1, py2, py3, py4]; % all y-coordinates

px = [2.02, 4.0, 6.00, 6.37];
py = [-0.41, -0.44, -0.589, 0.138];

umin = -0.45;
umax = -0.25;

umin = -0.40;
umax = -0.25;

figure
set(gcf,'Position',[100 100 width height])
contourf(x,y,field(:,:,1),200,'LineStyle','none')
hold on
scatter(px, py, 20,'filled','MarkerFaceColor',[255, 0, 0]/255,'MarkerFaceAlpha',0.7); % 'filled' just makes the markers solid
axis([xMinROI xMaxROI yMinROI yMaxROI])
colorbar
%colormap(mycolormap)
colormap(cmap)
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
fill(xcyl,ycyl,[1 1 1])  % place cylinder
plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
xlim([-2,20])
ylim([-4,4])
xticks(0:4:20)
yticks(-4:4:4)
%xlim([-1,8])
%ylim([-2,2])
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
drawnow

%%
load("bFTLE_original.mat")
%%
FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
ep=0.01;
SetFont = 'Times New Roman';
width =800; 
height = 320;
axislength = 20;
cmap = customcolormap(linspace(0,1,5),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000'}); %bluegreen
field = FTLE;
%field(FTLE<0.10) = nan;
% Example point coordinates
px = [px1, px2, px3, px4]; % all x-coordinates
py = [py1, py2, py3, py4]; % all y-coordinates

px = [2.02, 4.0, 6.00, 6.37];
py = [-0.41, -0.44, -0.589, 0.138];

umin = -0.45;
umax = -0.25;

umin = 0;
umax = 0.25;

figure
set(gcf,'Position',[100 100 width height])
contourf(x,y,field(:,:,1),200,'LineStyle','none')
hold on
%scatter(px, py, 20,'filled','MarkerFaceColor',[255, 0, 0]/255,'MarkerFaceAlpha',0.7); % 'filled' just makes the markers solid
axis([xMinROI xMaxROI yMinROI yMaxROI])
colorbar
%colormap(mycolormap)
colormap(cmap)
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
fill(xcyl,ycyl,[1 1 1])  % place cylinder
plot(xcyl,ycyl,'k','LineWidth',1.2) % cylinder boundary
xlim([-2,20])
ylim([-4,4])
xticks(0:4:20)
yticks(-4:4:4)
%xlim([-1,8])
%ylim([-2,2])
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
drawnow
