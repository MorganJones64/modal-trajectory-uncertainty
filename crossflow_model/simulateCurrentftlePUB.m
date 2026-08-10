%% Vector Field Trajectory Plus Particles Video
clear all
close all
addpath('models')
graphicsON = 1;   % flag for graphics
tstart = tic;
A = 0.1;    % parameters from Shadden 2005 Physica D
gamma = 0.25;
om = 0; %2*pi/5;  % frequency of gyre oscillations
omega = 2*pi/5;
bt=4;
at = 1.4;
eps = 0.03; %0.03; %0.1; %0.02;
dt =0.025;  % timestep
dx = 0.006; %0.015; %0.005; % %
dxv = 0.0125;
T = 15;     % duration of integration
int = 'f'; %''f for forward integration, 'b' for backward integration;
sf = 0.2; %3.4; %4;
xmax = 2; %4;  %vector range
ymax = 1; %2; %5;
xmin = 0;
ymin = -0; %-1;

xvec = xmin:dxv:xmax; %vector range
yvec = ymin:dxv:ymax;

xgmax = 2; %1.5; %Grid range
ygmax = 1; %1.5;
xgmin = 0;
ygmin = 0.;

xgrid = xgmin:dx:xgmax;
ygrid = ygmin:dx:ygmax;
[x0grid,y0grid] = meshgrid(xgrid,ygrid);  % grid of particles for later t
[x0,y0] = meshgrid(xgrid,ygrid);  % grid of particles for t=0
[x0vec,y0vec] = meshgrid(xvec,yvec);

yvIC(1,:,:) = x0vec'; %v is for vector display
yvIC(2,:,:) = y0vec';
yIC(1,:,:) = x0'; %for scatter point display
yIC(2,:,:) = y0';
if int == 'f'
    sgn = 1;
    tVec =0:(T/(dt));
else
    sgn = -1;
    tVec = flip(0:(T/(dt)));
end

color = 'k';
yin = yIC; %Show at the start (t=0) for video purposes
yout = yin; %Use to compute other points (t>0)


%set(gcf,'Position',[100 100 600 264.15])

step=8; %15;
color = linspace(1,10,length(xgrid)*length(ygrid));
%v = VideoWriter('currentclean_S9');
%v.FrameRate = 20;
%open(v)
%%
sf = 0.1; %0.2; %3.4; %4;
fontsize=15
if graphicsON
    figure(1)
    %set(gcf,'Position',[0 0 1000 1000])
    set(gcf,'Position',[0 0 1000 500])
end
for i=tVec
    time = i*dt;
        if graphicsON
            %subplot(2,1,1)
            %scatter(yout(1,:),yout(2,:),15,color,'o','filled','MarkerFaceAlpha',0.4)
            %colormap(gca,"jet")
            %axis off
            %axis equal
            %axis([xmin xmax ymin ymax])
            %set(gcf,'color','w')
            %yticks(linspace(-0.5,1.5,5))
            %ax=gca; ax.FontSize=fontsize;
            %colorbar0
            %subplot(2,1,2)
            dy = gVEC(time,yvIC,omega,bt,at);
            %dy = currentVEC(time,yvIC,A,gamma,om,0,omega,bt,at);
            quiver(yvIC(1,1:step:end,1:step:end),yvIC(2,1:step:end,1:step:end),sf*dy(1,1:step:end,1:step:end),sf*dy(2,1:step:end,1:step:end),'off','Color','#808080','linewidth',1.8);
            axis([xmin xmax ymin ymax])
            xticks(linspace(0,2,5))
            yticks(linspace(0,1,3))
            %xticklabels([])
            %yticklabels([])
            ax=gca; ax.FontSize=fontsize; %ax.FontName = 'Times New Roman'
            set(gcf,'color','w')
            box on
            drawnow
        end
        yout = rk4singlestep(@(t,x)currentVEC(t,x,A,gamma,om,eps,omega,bt,at),sgn*dt,sgn*time,yin);
        yin = yout;    
        %img = getframe(1);
        %writeVideo(v,img); 
        %hold on
        %hold off
end
%close(v)
% reshape 3-dim array into 2-dim array
xT = reshape(yout(1,:,:),length(xgrid),length(ygrid));
yT = reshape(yout(2,:,:),length(xgrid),length(ygrid));


%% Compute the finite-time Lyapunov exponent (sigma) Single Snapshot
% Finite difference to compute the g0radient
[dxTdx0,dxTdy0] = gradient(xT,dx,dx);
[dyTdx0,dyTdy0] = gradient(yT,dx,dx);
if int == 'f'
    mycolormap = customcolormap(linspace(0,1,7),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#ffffff'});
else
    mycolormap = customcolormap(linspace(0,1,6),{'#132d8d','#0d449d','#0758ab','#2295e0','#a9d9f8','#ffffff'});
end
% compute sigma: large sigma indicates large mixing!
for i=1:length(xgrid)
    for j=1:length(ygrid)
        D(1,1) = dxTdx0(i,j);
        D(1,2) = dxTdy0(i,j);
        D(2,1) = dyTdx0(i,j);
        D(2,2) = dyTdy0(i,j);
        sigma(i,j) = (1/T)*log(sqrt(max(eig(D'*D))));
    end
end

figure
set(gcf,'Position',[100 100 1000 500])
contourf(x0',y0',sigma,200,'LineStyle','none')
%set(gcf,'Position',[100 100 600 300])
axis([xgmin xgmax ygmin ygmax])
%colorbar
xticklabels([])
yticklabels([])
c=colorbar
c.FontSize =25;
c.Ticks = [0, 0.25, 0.5]
c.FontName = 'Times New Roman'
colormap(mycolormap)
clim([0,0.3])
set(gcf,'color','w')
%axis equal

%% Compute the finite-time Lyapunov exponent (sigma) Multi-Snapshot Video
xgrid = xgmin:dx:xgmax;
ygrid = ygmin:dx:ygmax;
[x0grid,y0grid] = meshgrid(xgrid,ygrid);  % grid of particles for later t

yIC(1,:,:) = x0grid';
yIC(2,:,:) = y0grid';

intLength = T/dt; %How long the integration is for
intDur = 250; %intLength/2; %How many FTLE flow fields should be saved
for r = 0:intDur
    r
    if int == 'f'
        sgn = 1;
        tVec = r:(T/dt)+(r);
    else
        sgn = -1;
        tVec = flip(r:(T/dt)+(r));
    end
    
    yin = yIC;
    for i=tVec
        time = i*dt;
        yout = rk4singlestep(@(t,x)currentVEC(t,x,A,gamma,om,0,omega,bt,at),sgn*dt,sgn*time,yin);
        yin = yout;   
    end
    % reshape 3-dim array into 2-dim array (final positions)
    xT = reshape(yout(1,:,:),length(xgrid),length(ygrid));
    yT = reshape(yout(2,:,:),length(xgrid),length(ygrid));

    [dxTdx0,dxTdy0] = gradient(xT,dx,dx);
    [dyTdx0,dyTdy0] = gradient(yT,dx,dx);

    % compute sigma: large sigma indicates large mixing!
    for i=1:length(xgrid)
        for j=1:length(ygrid)
            D(1,1) = dxTdx0(i,j);
            D(1,2) = dxTdy0(i,j);
            D(2,1) = dyTdx0(i,j);
            D(2,2) = dyTdy0(i,j);
            sigma(i,j,r+1) = max(eig(D'*D));
        end
    end
end

notes = 'eps=0.03 T=15s to match CSE parameters. NOTE: sigma is not scaled with sqrt or log or 1/T'
save("FTLEcurrent-smallamp_newperiod.mat")
