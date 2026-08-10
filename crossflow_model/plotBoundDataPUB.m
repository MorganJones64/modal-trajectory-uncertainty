%% Trajectory Errors p1
addpath('data')
%p1 = load('currentparticles_pert_location1_pert_b4_T24_point.mat');
%p2 = load('currentparticles_pert_location1_clean_b4_T24_point.mat');
%d1 = load('currentboundmse_tff_location1_b4_T24.mat');

%p1 = load('currentparticles_pert_location2_pert_b4_T24_point.mat');
%p2 = load('currentparticles_pert_location2_clean_b4_T24_point.mat');
%d1 = load('currentboundmse_tff_location2_b4_T24.mat');

%p1 = load('currentparticles_pert_location3_pert_b4_T24_point.mat');
%p2 = load('currentparticles_pert_location3_clean_b4_T24_point.mat');
%d1 = load('currentboundmse_tff_location3_b4_T24.mat');

%p1 = load('currentparticles_pert_location5_pert_b4_T24_point.mat');
%p2 = load('currentparticles_pert_location5_clean_b4_T24_point.mat');
%d1 = load('currentboundmse_tff_location5_b4_T24.mat');

p1 = load('currentparticles_pert_location6_pert_b4_T24_point.mat');
p2 = load('currentparticles_pert_location6_clean_b4_T24_point.mat');
d1 = load('currentboundmse_tff_location6_b4_T24.mat');
%%
%close all
ii = 3;
jj = 3;
SetFont = 'Times New Roman';
axislength = 18;
ep=0.01;
tLength = d1.intLength;
dt = d1.dt;

c1 = [44, 153, 61, 255]/255;
c2 = [89, 99, 235, 255]/255;
c3 = [235, 99, 89, 255]/255;
c4 = [203, 56, 232, 255]/255;

tVec = (0:1:tLength-1)*dt;
cseInt1 = reshape(d1.mseIntegral(ii,jj,:),[tLength,1]);
delta1 = reshape(d1.deltaInfty(ii,jj,:),[tLength,1]);
err1 = ((p1.xT-p2.xT).^2 + (p1.yT-p2.yT).^2).^2;

fig=figure('Position', [10 10 470 470])
h1 = semilogy(tVec,(cseInt1.*delta1).^2,'color',c2,'LineWidth',2,'LineStyle','-')
hold on
for i = 1:20:100
    for j = 1:20:100
        error1 = reshape(err1(j,i,:),[1,tLength]);
        semilogy(tVec,error1(1:tLength),'color',c4)      
    end
end
mean_error1 = reshape(mean(err1,[1,2]),[1,tLength]);
semilogy(tVec,mean_error1,'color','k','LineWidth',2) 
xlim([0,20])
%ylim([10^-8,10^1])
%ylim([10^-10,10^2])
ylim([10^-15,10^5])
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
xlabel('$t$','Interpreter','latex','fontsize',axislength);
ylabel('Squared Error','fontsize',axislength);
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
exportgraphics(fig, 'bound_location6.tif', 'Resolution', 1200)
%% Delta
figure('Position', [10 10 470 470])
h1 = semilogy(tVec,(delta1),'color',c2,'LineWidth',2,'LineStyle','-')
hold on
%h2 = semilogy(tVec,(cseInt1.*d1.epsn).^2,'color',c3,'LineWidth',2,'LineStyle','-')
for i = 1:20:100
    for j = 1:20:100
%        error1 = reshape(err1(j,i,:),[1,tLength]);
%        semilogy(tVec,error1(1:tLength),'color',c4)      
    end
end
mean_error1 = reshape(mean(err1,[1,2]),[1,tLength]);
%semilogy(tVec,mean_error1,'color','k','LineWidth',2) 
xlim([0,20])
%ylim([10^-8,10^1])
ylim([10^-10,10^2])
ylim([10^-20,10^2])
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
xlabel('$t$','Interpreter','latex','fontsize',axislength);
ylabel('Squared Error','fontsize',axislength);
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;

%% Plot all lines
% Trajectory Errors p1
addpath('data')
p1a = load('currentparticles_pert_location1_pert_b4_T24_point.mat');
p2a = load('currentparticles_pert_location1_clean_b4_T24_point.mat');
d1a = load('currentboundmse_tff_location1_b4_T24.mat');

p1b = load('currentparticles_pert_location2_pert_b4_T24_point.mat');
p2b = load('currentparticles_pert_location2_clean_b4_T24_point.mat');
d1b = load('currentboundmse_tff_location2_b4_T24.mat');

p1c = load('currentparticles_pert_location3_pert_b4_T24_point.mat');
p2c = load('currentparticles_pert_location3_clean_b4_T24_point.mat');
d1c = load('currentboundmse_tff_location3_b4_T24.mat');
%%
%close all
ii = 3;
jj = 3;
SetFont = 'Times New Roman';
axislength = 18;
ep=0.01;
tLength = size(p1.xT,3);
dt = d1.dt;

c1 = [44, 153, 61, 255]/255;
c2 = [89, 99, 235, 255]/255;
c3 = [235, 99, 89, 150]/255;
c4 = [203, 56, 232, 255]/255;

tVec = (0:1:tLength-1)*dt;
cseInta = reshape(d1a.mseIntegral(ii,jj,:),[tLength,1]);
deltaa = reshape(d1a.deltaInfty(ii,jj,:),[tLength,1]);
erra = ((p1a.xT-p2a.xT).^2 + (p1a.yT-p2a.yT).^2).^2;

cseIntb = reshape(d1b.mseIntegral(ii,jj,:),[tLength,1]);
deltab = reshape(d1b.deltaInfty(ii,jj,:),[tLength,1]);
errb = ((p1b.xT-p2b.xT).^2 + (p1b.yT-p2b.yT).^2).^2;

cseIntc = reshape(d1c.mseIntegral(ii,jj,:),[tLength,1]);
deltac = reshape(d1c.deltaInfty(ii,jj,:),[tLength,1]);
errc = ((p1c.xT-p2c.xT).^2 + (p1c.yT-p2c.yT).^2).^2;

figure('Position', [10 10 470 470])
ha = semilogy(tVec,(cseInta.*deltaa).^2,'color',c1,'LineWidth',2,'LineStyle','-')
hold on
hb = semilogy(tVec,(cseIntb.*deltab).^2,'color',c2,'LineWidth',2,'LineStyle','-')
%hc = semilogy(tVec,(cseIntc.*deltac).^2,'color',c3,'LineWidth',3,'LineStyle','-')
for i = 1:20:100
    for j = 1:20:100
        errora = reshape(erra(j,i,:),[1,tLength]);
        semilogy(tVec,errora(1:tLength),'color',c1,'LineStyle','-.','LineWidth',2)  
        errorb = reshape(errb(j,i,:),[1,tLength]);
        semilogy(tVec,errorb(1:tLength),'color',c2,'LineStyle','-.','LineWidth',2)  
        %errorc = reshape(errc(j,i,:),[1,tLength]);
        %semilogy(tVec,errorc(1:tLength),'color',c3)  
    end
end
xlim([0,20])
%ylim([10^-8,10^1])
%ylim([10^-10,10^2])
ylim([10^-15,10^5])
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
xlabel('$t$','Interpreter','latex','fontsize',axislength);
ylabel('Squared Error','fontsize',axislength);
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;
