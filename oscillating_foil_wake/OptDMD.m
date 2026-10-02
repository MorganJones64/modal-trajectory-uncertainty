clear all
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
addpath('functions\optdmdsrc\')
addpath('customcolormap')
load('data\oscfoil_data.mat')
 %% Preprocessing Flow field data `(Temporal DMD)
U(:,:,:) = u;
U(:,:,nt+1:2*nt) = v;
%%
% Reshape the velocity matrix for DMD
U = [reshape(u,[nx*ny,nt]); ...
     reshape(v,[nx*ny,nt])];
%%
r=9;
imode=1;
%Apply constraints to only allow a periodic mode. (No growth/decay)
      %real part of alpha ;  imaginary part of alpha in the rectangular
      %space
lbc = [zeros([r, 1]); -Inf*ones([r, 1])];
ubc = [zeros([r, 1]); Inf*ones([r, 1])];


copts = varpro_lsqlinopts('lbc',lbc,'ubc',ubc);
%%
[w,e,b] = optdmd(U,tVec,r,imode,[],[],[],copts);
%% sort mode coefficients based on frequency
[~,ind]=sort(abs(imag(e)),'ascend')
b=b(ind);
e=e(ind);
w=w(:,ind);
gamma = exp(e*dt);
save('data\optdmd_coeffs.mat', 'w','e','b','gamma','r','copts')

%% Plot amplitude coeffients Figure 9a
fontSize = 20;
ep=0.01
width =450;
height = 350;
A=0.02227;
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
set(gca,'LineWidth',1.2,'TickLength',[ep, 0],'YMinorGrid','off','XMinorGrid','off')
box on
set(gcf,'color','white')
grid off

%% Plot Eigenvalues Figure 9b
close all
ep=0.01
fontSize=20
width =450;
height = 350;
%r=indsLambda(1);
figure(10)
set(10,'position',[200,0,width,height])
setFont = 'Times New Roman';
hold on
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k--')

plot(real(gamma(1)),imag(gamma(1)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',lw)
plot(real(gamma(2:3)),imag(gamma(2:3)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c4,'LineWidth',lw)
plot(real(gamma(4:5)),imag(gamma(4:5)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',c1,'LineWidth',lw)
plot(real(gamma(6:7)),imag(gamma(6:7)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',green,'LineWidth',lw)
plot(real(gamma(8:9)),imag(gamma(8:9)),'o','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',blue,'LineWidth',lw)
hold on

xlabel('$\Re(\lambda_i)$','interpreter','latex','fontsize',20)
ylabel('$\Im(\lambda_i)$','interpreter','latex','fontsize',20)
ax=gca; ax.FontName = setFont; ax.FontSize = fontSize; 
box on
axis equal
set(gca,'LineWidth',1.2,'TickLength',[ep, 0],'YMinorGrid','off','XMinorGrid','off')
set(gcf,'color','white')
xlim([0.8 1.2])
ylim([-0.40 0.40])
yticks([linspace(-0.4,0.4,3)])
grid off

%% Compute modal structures
Umean = real(w(:,1)*diag(b(1))*exp(e(1)*tVec));
U1 = real(w(:,2:3)*diag(b(2:3))*exp(e(2:3)*tVec));
U2 = real(w(:,4:5)*diag(b(4:5))*exp(e(4:5)*tVec));
U3 = real(w(:,6:7)*diag(b(6:7))*exp(e(6:7)*tVec));
U4 = real(w(:,8:9)*diag(b(8:9))*exp(e(8:9)*tVec));

nu = 1:nx*ny;
nv = (nx*ny)+1:2*(nx*ny);

umean = reshape(Umean(nu,:),[ny nx nt]);
vmean = reshape(Umean(nv,:),[ny nx nt]);

u1 = reshape(U1(nu,:),[ny nx nt]);
v1 = reshape(U1(nv,:),[ny nx nt]);

u2 = reshape(U2(nu,:),[ny nx nt]);
v2 = reshape(U2(nv,:),[ny nx nt]);

u3 = reshape(U3(nu,:),[ny nx nt]);
v3 = reshape(U3(nv,:),[ny nx nt]);

u4 = reshape(U4(nu,:),[ny nx nt]);
v4 = reshape(U4(nv,:),[ny nx nt]);

save('data/opdmd_modes.mat', 'umean','u1','u2','u3','u4','vmean','v1','v2','v3','v4')

%% Dyn Sys Full Size Original Figure 8
close all
addpath('./customcolormap/')
load('foilcoords.mat')
mycolormap = customcolormap(linspace(0,1,13),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});
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
ind = [4];
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
i=period_sq([4])%:800
% 
fig1=figure
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


fig2=figure
umin = 0.5;
umax = 1.5;
set(gcf,'Position',[100 100 width height])
uplot = u1(:,:,i) + umean(:,:,i);
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



