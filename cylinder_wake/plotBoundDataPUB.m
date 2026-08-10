px1 = 2.013;
py1 = -0.424;

px2 = 4.183;
py2 = -0.424;

px3 = 4.806;
py3 = 0.869;

px4 = 6.0;
py4 = 0.061;


%1-4 in order of increasing x location
%% Trajectory Errors m1-2
%close all
clear all
floc = './data/';
addpath('./data')
p41 = load([floc,'particlePos_backwards_p8_original.mat']);
p42 = load([floc,'particlePos_backwards_p8_fmnm1m2.mat']);

px4 = p41.xMinROI + p41.dtol;
py4 = p41.yMinROI + p41.dtol; 

p31 = load([floc,'particlePos_backwards_p7_original.mat']);
p32 = load([floc,'particlePos_backwards_p7_fmnm1m2.mat']);

px3 = p31.xMinROI + p31.dtol;
py3 = p31.yMinROI + p31.dtol; 

p21 = load([floc,'particlePos_backwards_p6b_original.mat']);
p22 = load([floc,'particlePos_backwards_p6b_fmnm1m2.mat']);

px2 = p21.xMinROI + p21.dtol;
py2 = p21.yMinROI + p21.dtol; 

p11 = load([floc,'particlePos_backwards_p5_original.mat']);
p12 = load([floc,'particlePos_backwards_p5_fmnm1m2.mat']);

px1 = p11.xMinROI + p11.dtol;
py1 = p11.yMinROI + p11.dtol; 

d4 = load('bounddatasave_backwards_p8_fmnm1m2_gm3.mat');
d3 = load('bounddatasave_backwards_p7_fmnm1m2_gm3.mat');
d2 = load('bounddatasave_backwards_p6b_fmnm1m2_gm3.mat');
d1 = load('bounddatasave_backwards_p5_fmnm1m2_gm3.mat');

%
%close all
ii = 2;
jj = 2;
SetFont = 'Times New Roman';
axislength = 18;
ep=0.01;
tLength = d4.tLength;
dt = d4.dt

%% mnm1m2
%close all
%c2 = [89, 99, 235, 120]/255;
c2o = [89, 99, 235, 255]/255;
%c1 = [44, 153, 61, 120]/255;
c1o = [44, 153, 61, 255]/255;
%c3 = [235, 99, 89, 120]/255;
c3o = [235, 99, 89, 255]/255;
%c4 = [235, 99, 89, 120]/255;
%c4o = [235, 99, 89, 255]/255;

tVec = (0:1:abs(tLength)-1)/60;
cseInt4 = reshape(d4.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta4 = reshape(d4.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta4 = reshape((1/abs(tLength*dt))*(log(d4.deltaInfty(ii,jj,:).*d4.cseIntegral(ii,jj,:)./p42.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

cseInt3 = reshape(d3.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta3 = reshape(d3.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta3 = reshape((1/abs(tLength*dt))*(log(d3.deltaInfty(ii,jj,:).*d3.cseIntegral(ii,jj,:)./p32.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

cseInt2 = reshape(d2.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta2 = reshape(d2.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta2 = reshape((1/abs(tLength*dt))*(log(d2.deltaInfty(ii,jj,:).*d2.cseIntegral(ii,jj,:)./p22.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

cseInt1 = reshape(d1.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta1 = reshape(d1.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta1 = reshape((1/abs(tLength*dt))*(log(d1.deltaInfty(ii,jj,:).*d1.cseIntegral(ii,jj,:)./p12.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

err4 = (p41.xPos-p42.xPos).^2 + (p41.yPos-p42.yPos).^2;
err3 = (p31.xPos-p32.xPos).^2 + (p31.yPos-p32.yPos).^2;
err2 = (p21.xPos-p22.xPos).^2 + (p21.yPos-p22.yPos).^2;
err1 = (p11.xPos-p12.xPos).^2 + (p11.yPos-p12.yPos).^2;

%figure('Position', [10 10 470 470])
figure('Position', [10 10 470 500])
h4 = semilogy(tVec,(cseInt4.*delta4),'color','k','LineWidth',2,'LineStyle','-');
hold on
h3 = semilogy(tVec,(cseInt3.*delta3),'color',c1o,'LineWidth',2,'LineStyle','-');
h2 = semilogy(tVec,(cseInt2.*delta2),'color',c2o,'LineWidth',2,'LineStyle','-');
h1 = semilogy(tVec,(cseInt1.*delta1),'color',c3o,'LineWidth',2,'LineStyle','-');
hold on
for i = 2
    for j = 2
        semilogy(tVec,reshape(err4(i,j,1:abs(tLength)),[1 abs(tLength)]),'color','k','Linewidth',2,'LineStyle',':')
        semilogy(tVec,reshape(err3(i,j,1:abs(tLength)),[1 abs(tLength)]),'color',c1o,'Linewidth',2,'LineStyle',':')
        semilogy(tVec,reshape(err2(i,j,1:abs(tLength)),[1 abs(tLength)]),'color',c2o,'Linewidth',2,'LineStyle',':')
        semilogy(tVec,reshape(err1(i,j,1:abs(tLength)),[1 abs(tLength)]),'color',c3o,'Linewidth',2,'LineStyle',':')
    end
end

xlim([0,1])
ylim([10^-6,10^0])
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
xlabel('$t/T_{p}$','Interpreter','latex','fontsize',axislength);
ylabel('Squared Error','fontsize',axislength);

legend([h4 h3 h2 h1], ...
{
sprintf('$(x,y) = (%.2f, %.2f)$', px4, py4), ...
sprintf('$(x,y) = (%.2f, %.2f)$', px3, py3), ...
sprintf('$(x,y) = (%.2f, %.2f)$', px2, py2), ...
sprintf('$(x,y) = (%.2f, %.2f)$', px1, py1)
}, ...
'Location','southeast', ...
'FontSize',15, ...
'Interpreter','latex')

%legend([h1 h2 h3 h4], {'region-1','region-2','region-3','region-4'}, ...
%       'location','southeast')
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;

%% Zeta
figure('Position', [10 10 470 470])
h4 = plot(tVec,zeta4,'color','k','LineWidth',2,'LineStyle','-');
hold on
h3 = plot(tVec,zeta3,'color',c1o,'LineWidth',2,'LineStyle','-');
h2 = plot(tVec,zeta2,'color',c2o,'LineWidth',2,'LineStyle','-');
h1 = plot(tVec,zeta1,'color',c3o,'LineWidth',2,'LineStyle','-');

set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
xlabel('$t/T_{p}$','Interpreter','latex','fontsize',axislength);
ylabel('\zeta','fontsize',axislength);
legend([h4 h3 h2 h1], ...
       {'$p8-(x,y) = (6.37, 0.14)$', ...
        '$p7-(x,y) = (6.01, -0.59)$', ...
        '$p6-(x,y) = (4.17, -0.44)$', ...
        '$p5-(x,y) = (2.02, -0.41)$'}, ...
       'Location','southeast','FontSize',15,'Interpreter','latex')
%%
[zeta1(1),zeta2(1),zeta3(1),zeta4(1)]

[zeta1(end-2),zeta2(end-2),zeta3(end-2),zeta4(end-2)]


%% Farfield region Trajectory Errors
%close all
clear all
floc = './data/boundextra/';
floc2 = './data/';
addpath('./data')
addpath('./data/boundextra/')
p41 = load([floc,'particlePos_backwards_p11_original.mat']);
p42 = load([floc,'particlePos_backwards_p11_fmnm1m2.mat']);

px4 = p41.xMinROI + p41.dtol
py4 = p41.yMinROI + p41.dtol 

%p10
p31 = load([floc,'particlePos_backwards_p13_original.mat']);
p32 = load([floc,'particlePos_backwards_p13_fmnm1m2.mat']);
px3 = p31.xMinROI + p31.dtol
py3 = p31.yMinROI + p31.dtol

p21 = load([floc2,'particlePos_backwards_p6c_original.mat']);
p22 = load([floc2,'particlePos_backwards_p6c_fmnm1m2.mat']);

px2 = p21.xMinROI + p21.dtol
py2 = p21.yMinROI + p21.dtol

p11 = load([floc2,'particlePos_backwards_p5_original.mat']);
p12 = load([floc2,'particlePos_backwards_p5_fmnm1m2.mat']);

px1 = p11.xMinROI + p11.dtol
py1 = p11.yMinROI + p11.dtol 

d4 = load('bounddatasave_backwards_p11_fmnm1m2_gm3.mat');
d3 = load('bounddatasave_backwards_p13_fmnm1m2_gm3.mat');
d2 = load('bounddatasave_backwards_p6c_fmnm1m2_gm3.mat');
d1 = load('bounddatasave_backwards_p5_fmnm1m2_gm3.mat');

%
%close all
ii = 2;
jj = 2;
SetFont = 'Times New Roman';
axislength = 18;
ep=0.01;
tLength = d4.tLength;
dt = d4.dt

%% mnm1m2
%close all
%c2 = [89, 99, 235, 120]/255;
c2o = [89, 99, 235, 255]/255;
%c1 = [44, 153, 61, 120]/255;
c1o = [44, 153, 61, 255]/255;
%c3 = [235, 99, 89, 120]/255;
c3o = [235, 99, 89, 255]/255;
%c4 = [235, 99, 89, 120]/255;
%c4o = [235, 99, 89, 255]/255;

tVec = (0:1:abs(tLength)-1)/60;
cseInt4 = reshape(d4.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta4 = reshape(d4.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta4 = reshape((1/abs(tLength*dt))*(log(d4.deltaInfty(ii,jj,:).*d4.cseIntegral(ii,jj,:)./p42.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

cseInt3 = reshape(d3.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta3 = reshape(d3.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta3 = reshape((1/abs(tLength*dt))*(log(d3.deltaInfty(ii,jj,:).*d3.cseIntegral(ii,jj,:)./p32.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

cseInt2 = reshape(d2.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta2 = reshape(d2.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta2 = reshape((1/abs(tLength*dt))*(log(d2.deltaInfty(ii,jj,:).*d2.cseIntegral(ii,jj,:)./p22.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

cseInt1 = reshape(d1.cseIntegral(ii,jj,:),[abs(tLength),1]);
delta1 = reshape(d1.deltaInfty(ii,jj,:),[abs(tLength),1]);
zeta1 = reshape((1/abs(tLength*dt))*(log(d1.deltaInfty(ii,jj,:).*d1.cseIntegral(ii,jj,:)./p12.sigma_ftle(ii,jj,:))),[abs(tLength),1]);;

err4 = (p41.xPos-p42.xPos).^2 + (p41.yPos-p42.yPos).^2;
err3 = (p31.xPos-p32.xPos).^2 + (p31.yPos-p32.yPos).^2;
err2 = (p21.xPos-p22.xPos).^2 + (p21.yPos-p22.yPos).^2;
err1 = (p11.xPos-p12.xPos).^2 + (p11.yPos-p12.yPos).^2;

%figure('Position', [10 10 470 470])
figure('Position', [10 10 470 500])
h4 = semilogy(tVec,(cseInt4.*delta4),'color','k','LineWidth',2,'LineStyle','-');
hold on
h3 = semilogy(tVec,(cseInt3.*delta3),'color',c2o,'LineWidth',2,'LineStyle','-');
h2 = semilogy(tVec,(cseInt2.*delta2),'color',c3o,'LineWidth',2,'LineStyle','-');
%h1 = semilogy(tVec,(cseInt1.*delta1),'color',c1o,'LineWidth',2,'LineStyle','-');
hold on
for i = 2
    for j = 2
        semilogy(tVec,reshape(err4(i,j,1:abs(tLength)),[1 abs(tLength)]),'color','k','Linewidth',2,'LineStyle',':')
        semilogy(tVec,reshape(err3(i,j,1:abs(tLength)),[1 abs(tLength)]),'color',c2o,'Linewidth',2,'LineStyle',':')
        semilogy(tVec,reshape(err2(i,j,1:abs(tLength)),[1 abs(tLength)]),'color',c3o,'Linewidth',2,'LineStyle',':')
        %semilogy(tVec,reshape(err1(i,j,1:abs(tLength)),[1 abs(tLength)]),'color',c1o,'Linewidth',2,'LineStyle',':')
    end
end

xlim([0,1])
%ylim([10^-7,10^1])
ylim([10^-8,10^0])
%yticks([10^-7 10^-5 10^-3 10^-1 10^1])
set(gca,'LineWidth',1.2,'TickLength',[ep, ep])
set(gcf,'color','w')
xlabel('$t/T_{p}$','Interpreter','latex','fontsize',axislength);
ylabel('Squared Error','fontsize',axislength);

legend([h4 h3 h2], ...
{
sprintf('$(x,y) = (%.2f, %.2f)$', px4, py4), ...
sprintf('$(x,y) = (%.2f, %.2f)$', px3, py3), ...
sprintf('$(x,y) = (%.2f, %.2f)$', px2, py2), ...
%sprintf('$(x,y) = (%.2f, %.2f)$', px1, py1)
}, ...
'Location','southeast', ...
'FontSize',18, ...
'Interpreter','latex')

%legend([h1 h2 h3 h4], {'region-1','region-2','region-3','region-4'}, ...
%       'location','southeast')
ax2 = gca; box on; ax2.FontName = SetFont; ax2.FontSize = axislength;

%%
load("CSE_fm1_2_gm3_DNS_800_400.mat")
%%
xq = [px4, px3, px2 px1];
yq = [py4 , py3, py2 py1];
zeta = (1/abs(tLength*dt))*(log(deltaInfty)+log(cseIntegral)-log(sigma_ftle));
FTLE = (1/abs(tLength*dt))*log(sigma_ftle);
xMat = xPos(:,:,1);
yMat = yPos(:,:,1);
F = scatteredInterpolant( ...
    xMat(:), ...
    yMat(:), ...
    zeta(:), ...
    'linear', ...
    'nearest');

vq = F(xq,yq)
