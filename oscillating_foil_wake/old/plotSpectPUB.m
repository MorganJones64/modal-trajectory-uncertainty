%save('DNS_DMDtempcirc_int_spectcn.mat','ainds','amp','lambda','omega','dx','dt','dy','num','note')

clear all
addpath('data\')
load("optdmd_coeff.mat")

%% Compute reconstructions using modes with largest |lambda|
fontSize = 20;
ep=0.01
width =450;
height = 350;
A=0.02227;
St = A*imag(e)/(0.1*pi);
orange = [1, 0.6039, 0.2784];
yellow = [1, 0.9412, 0.2784];
blue = [0.1059, 0.4000, 0.9804];
green = [0, 0.7020, 0.2118];
c1='#f0b73d'; %'#f0993d';
c2='#eb4546';
c3='#da3f2b';
c4='#cf6548'; %'#da3f2b';
lw = 1.25;
lw2=1;
figure
set(gcf,'position',[200,0,width,height])
plot(imag(e(1)),b(1),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',lw)
hold on
plot(imag(e(2:3)),b(2:3),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c4,'LineWidth',lw)
plot(imag(e(4:5)),b(4:5),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c1,'LineWidth',lw)
plot(imag(e(6:7)),b(6:7),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',green,'LineWidth',lw)
plot(imag(e(8:9)),b(8:9),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',blue,'LineWidth',lw)
xlabel('$\Im(\omega_k$)','interpreter','latex','fontsize',fontSize)
ylabel('$b_k$','interpreter','latex','fontsize',12)
ax=gca; ax.FontName = 'Times New Roman'; ax.FontSize = fontSize; 
%xticks([linspace(-1.5,1.5,3)])
set(gca,'LineWidth',1.2,'TickLength',[ep, 0],'YMinorGrid','off','XMinorGrid','off')
box on
set(gcf,'color','white')
grid off



%% Plot Eigenvalues
close all
size = 70; %25;
ep=0.01
fontSize=20
width =450;
height = 350;
%r=indsLambda(1);
figure(10)
set(10,'position',[200,0,width,height])
setFont = 'Times New Roman';
%scatter(real(lambda),imag(lambda),'MarkerEdgeColor','b',...
%    'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0)
%alpha(s,0.5)
hold on
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k--')

plot(real(e(1)),imag(e(1)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',lw)
plot(real(e(2:3)),imag(e(2:3)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c4,'LineWidth',lw)
plot(real(e(4:5)),imag(e(4:5)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c1,'LineWidth',lw)
plot(real(e(6:7)),imag(e(6:7)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',green,'LineWidth',lw)
plot(real(e(8:9)),imag(e(8:9)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',blue,'LineWidth',lw)
hold on

%ylim([-0.8 0.8])
xlabel('$\Re(\lambda_i)$','interpreter','latex','fontsize',20)
ylabel('$\Im(\lambda_i)$','interpreter','latex','fontsize',20)
ax=gca; ax.FontName = setFont; ax.FontSize = fontSize; 
box on
axis equal
set(gca,'LineWidth',1.2,'TickLength',[ep, 0],'YMinorGrid','off','XMinorGrid','off')
set(gcf,'color','white')
%xlim([0 1.2])
xlim([0.8 1.2])
ylim([-0.40 0.40])
yticks([linspace(-0.4,0.4,3)])
%yticks(linspace(-0.15,0.15,3))
%xticks([0.9, 1, 1.1])
%pbaspect([1 1 1])
grid off
%title(['Temporal DMD Spectrum, St=',num2str(St)])

%% DARK  Compute reconstructions using modes with largest |lambda|
fontSize = 26;
ep=0.03
width =400;
height = 320;
A=0.02227;
St = A*imag(e)/(0.1*pi);
orange = [1, 0.6039, 0.2784];
yellow = [1, 0.9412, 0.2784];
blue = [0.1059, 0.4000, 0.9804];
green = [0, 0.7020, 0.2118];
c1='#f0b73d'; %'#f0993d';
c2='#eb4546';
c3='#da3f2b';
c4='#cf6548'; %'#da3f2b';
lw = 1.25;
lw2=1;
figure
set(gcf,'position',[200,0,width,height])
plot(imag(e(1)),b(1),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',lw)
hold on
plot(imag(e(2:3)),b(2:3),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',lw)
plot(imag(e(4:5)),b(4:5),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',lw)
plot(imag(e(6:7)),b(6:7),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',lw)
plot(imag(e(8:9)),b(8:9),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','w','LineWidth',lw)
% plot(St(1),b(1),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c4,'LineWidth',lw)
% hold on
% plot(St(2:3),b(2:3),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',orange,'LineWidth',lw)
% plot(St(4:5),b(4:5),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',yellow,'LineWidth',lw)
% plot(St(6:7),b(6:7),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',green,'LineWidth',lw)
% plot(St(8:9),b(8:9),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',blue,'LineWidth',lw)
%xlim([0,max(omega)])
%xlim([-1.5,1.5])
%ylim([0,250])
%xlabel('$\Im(St_k$)','interpreter','latex','fontsize',fontSize)
xlabel('$\Im(\omega_k$)','interpreter','latex','fontsize',12)
ylabel('$||a_k||_2$','interpreter','latex','fontsize',12)
ax=gca; ax.FontName = 'Times New Roman'; ax.FontSize = fontSize; ax.Color = 'k'
set(gca,'LineWidth',1.2,'TickLength',[ep, 0],'YMinorGrid','off','XMinorGrid','off')
ax.XColor = 'w'
ax.YColor = 'w'
box on
set(gcf,'color','k')
grid off

