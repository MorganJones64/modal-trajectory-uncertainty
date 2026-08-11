clear all
addpath('data\')
addpath('fxn\')
load('cylinderPOD.mat')
%% Singular Values and Cumulative Sum
close all
ep=0.01
sv1 = [1,2];
sv2 = [3,4];
sv3 = [5,6];
sv4 = [7,8];

width =800;
height = 450;
lw = 1.25;
lw2=1.8;

orange = [1, 0.6039, 0.2784];
yellow = [1, 0.9412, 0.2784];
blue = [0.1059, 0.4000, 0.9804];
green = [0, 0.7020, 0.2118];
c1='#f0b73d'; %'#f0993d';
c2='#eb4546';
c3='#da3f2b';
c4='#cf6548'; %'#da3f2b';

d = 40;
axislength=20;
SetFont = 'Times New Roman';
r=9;
figure
set(gcf,'position',[200,0,width,height])
subplot(1,2,1)
semilogy(sigma,'ok','LineWidth',lw)
hold on
semilogy(sv1,sigma(sv1),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c4,'LineWidth',lw)
semilogy(sv2,sigma(sv2),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c1,'LineWidth',lw)
semilogy(sv3,sigma(sv3),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',green,'LineWidth',lw)
semilogy(sv4,sigma(sv4),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',blue,'LineWidth',lw)
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
ylabel('$\sigma_{k}$','interpreter','latex')
xlabel('$k$','interpreter','latex')
%xlim([0,40])
set(gca,'LineWidth',1.2,'TickLength',[ep, 0],'YMinorGrid','off','XMinorGrid','off')
%grid on
hold off
grid on

csum = cumsum(sigma)/sum(sigma);
subplot(1,2,2) 
plot(cumsum(sigma)/sum(sigma),'ok','LineWidth',lw)
hold on
sc1 = cumsum(sigma(1:sv1))/sum(sigma);
sc2 = cumsum(sigma(1:sv2))/sum(sigma);
sc3 = cumsum(sigma(1:sv3))/sum(sigma);
plot(sv1,csum(sv1),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c4,'LineWidth',lw)
plot(sv2,csum(sv2),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c1,'LineWidth',lw)
plot(sv3,csum(sv3),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',green,'LineWidth',lw)
plot(sv4,csum(sv4),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',blue,'LineWidth',lw)

ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
%ylabel('\bf{Cumulative energy}','interpreter','latex')
%ylabel('$\sum_{n=1}^k\sigma_n/\sum_{n=1}^{n_t}\sigma_n$','interpreter','latex')
xlabel('$k$','interpreter','latex')
%xlim([0,40])
ylim([0.3,1])
yticks([0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
set(gcf,'Color','w')
set(gca,'LineWidth',1.2,'TickLength',[ep, 0])
grid on
