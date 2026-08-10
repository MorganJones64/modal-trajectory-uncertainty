clear all
addpath('data\')
addpath('fxn\')
load('cylindervelocity_zen.mat')
%%
umean = mean(u,3);
vmean = mean(v,3);
up = u - umean;
vp = v - vmean;
[Zeta, Sigma, Xi] = POD(up, vp, dx, dy, dt);

%% Singular Values and Cumulative Sum
sv1 = 2;
sv2 = 4;
sv3 = 8;
d = 40;
r=9;
figure
sigma = diag(Sigma);
subplot(1,2,1)
semilogy(sigma,'ok')
hold on
semilogy(sv1,sigma(sv1),'.r','MarkerSize',25)
text(sv1+d,sigma(sv1),['r =',num2str(sv1)],'interpreter','latex','fontsize',14,'Color','r')
semilogy(sv2,sigma(sv2),'.r','MarkerSize',25)
text(sv2+d,sigma(sv2),['r =',num2str(sv2)],'interpreter','latex','fontsize',14,'Color','r')
semilogy(sv3,sigma(sv3),'.r','MarkerSize',25)
text(sv3+d,sigma(sv3),['r =',num2str(sv3)],'interpreter','latex','fontsize',14,'Color','r')
ylabel('\bf{Singular value}, $\sigma{k}$','interpreter','latex')
xlabel('k','interpreter','latex')
%ylim([1e-5,1])
grid on
hold off

subplot(1,2,2) 
plot(cumsum(sigma)/sum(sigma),'ok')
hold on
sc1 = cumsum(sigma(1:sv1))/sum(sigma);
sc2 = cumsum(sigma(1:sv2))/sum(sigma);
sc3 = cumsum(sigma(1:sv3))/sum(sigma);
plot(sv1,sc1(end),'.r','MarkerSize',25)
text(sv1+d,sc1(end),['r =',num2str(sv1)],'interpreter','latex','fontsize',14,'Color','r')
plot(sv2,sc2(end),'.r','MarkerSize',25)
text(sv2+d,sc2(end),['r =',num2str(sv2)],'interpreter','latex','fontsize',14,'Color','r')
plot(sv3,sc3(end),'.r','MarkerSize',25)
text(sv3+d,sc3(end),['r =',num2str(sv3)],'interpreter','latex','fontsize',14,'Color','r')
%plot(cumsum(sigma(1:sv2))/sum(sigma(1:sv2)),'.r','MarkerSize',25)
%plot(cumsum(sigma(1:sv3))/sum(sigma(1:sv3)),'.r','MarkerSize',25)
ylabel('\bf{Cumulative energy}','interpreter','latex')
xlabel('k','interpreter','latex')
set(gcf,'Color','w')
grid on

%%

nu = 1:(nx*ny);
nv = (nx*ny)+1:2*(nx*ny);
num = 9:nt;
U1=Zeta(:,1:2)*Sigma(1:2,1:2)*Xi(:,1:2)';
U2=Zeta(:,3:4)*Sigma(3:4,3:4)*Xi(:,3:4)';
U3=Zeta(:,5:6)*Sigma(5:6,5:6)*Xi(:,5:6)';
U4=Zeta(:,7:8)*Sigma(7:8,7:8)*Xi(:,7:8)';
Ut=Zeta(:,num)*Sigma(num,num)*Xi(:,num)';
Umean = Zeta(:,1)*Sigma(1,1)*Xi(:,1)';

U11=Zeta(:,1)*Sigma(1,1)*Xi(:,1)';
U12=Zeta(:,2)*Sigma(2,2)*Xi(:,2)';

U21=Zeta(:,3)*Sigma(3,3)*Xi(:,3)';
U22=Zeta(:,4)*Sigma(4,4)*Xi(:,4)';

U31=Zeta(:,5)*Sigma(5,5)*Xi(:,5)';
U32=Zeta(:,6)*Sigma(6,6)*Xi(:,6)';

U41=Zeta(:,7)*Sigma(7,7)*Xi(:,7)';
U42=Zeta(:,8)*Sigma(8,8)*Xi(:,8)';

u11 = reshape(U11(nu,:),[ny nx nt]);
v11 = reshape(U11(nv,:),[ny nx nt]);
u12 = reshape(U12(nu,:),[ny nx nt]);
v12 = reshape(U12(nv,:),[ny nx nt]);

u21 = reshape(U21(nu,:),[ny nx nt]);
v21 = reshape(U21(nv,:),[ny nx nt]);
u22 = reshape(U22(nu,:),[ny nx nt]);
v22 = reshape(U22(nv,:),[ny nx nt]);

u31 = reshape(U31(nu,:),[ny nx nt]);
v31 = reshape(U31(nv,:),[ny nx nt]);
u32 = reshape(U32(nu,:),[ny nx nt]);
v32 = reshape(U32(nv,:),[ny nx nt]);

u41 = reshape(U41(nu,:),[ny nx nt]);
v41 = reshape(U41(nv,:),[ny nx nt]);
u42 = reshape(U42(nu,:),[ny nx nt]);
v42 = reshape(U42(nv,:),[ny nx nt]);

%ut = reshape(Ut(nu,:),[ny nx nt]);
%vt = reshape(Ut(nv,:),[ny nx nt]);

save('cylinderPOD.mat','Zeta','Sigma','Xi')
save('cylinderPODmodesSingle.mat','u11','u21','u31','u41','v11','v21','v31','v41','u12','u22','u32','u42','v12','v22','v32','v42','umean','vmean')

clear u11 u21 u31 u41 v11 v21 v31 v41 u12 u22 u32 u42 v12 v22 v32 v42

u1 = reshape(U1(nu,:),[ny nx nt]);
u2 = reshape(U2(nu,:),[ny nx nt]);
u3 = reshape(U3(nu,:),[ny nx nt]);
u4 = reshape(U4(nu,:),[ny nx nt]);

v1 = reshape(U1(nv,:),[ny nx nt]);
v2 = reshape(U2(nv,:),[ny nx nt]);
v3 = reshape(U3(nv,:),[ny nx nt]);
v4 = reshape(U4(nv,:),[ny nx nt]);

ut = reshape(Ut(nu,:),[ny nx nt]);
vt = reshape(Ut(nv,:),[ny nx nt]);
save('cylinderPODmodesPair.mat','u1','u2','u3','u4','ut','v1','v2','v3','v4','vt','umean','vmean')