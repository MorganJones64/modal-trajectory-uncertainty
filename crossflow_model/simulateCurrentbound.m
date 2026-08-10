clear all
close all
addpath('models')

A = 0.1;    % parameters from Shadden 2005 Physica D
gamma = 0.25;
U=0.2;
om = 0; %2*pi/5;  % frequency of gyre oscillations
omega = 2*pi/5;
bt=4;
at = 1.4;
epsn = 0.03; %0.1; %0.02;
dt =0.025;  % timestep
T = 15;     % duration of integration
int = 'f'; %''f for forward integration, 'b' for backward integration;

%location 1
% tol = 1e-7;
% xgmax = 0.48+tol; %Grid range
% ygmax = 0.63+tol; 
% xgmin = 0.48-tol;
% ygmin = 0.63-tol;

%location 2
tol = 1e-7;
xgmax = 1.5+tol; %Grid range
ygmax = 0.75+tol;
xgmin = 1.5-tol;
ygmin = 0.75-tol;

%location 3
% tol = 1e-7;
% xgmax = 1.0265; % Particle Grid range
% ygmax = 0.0665; 
% xgmin = 1.026;
% ygmin = 0.066;

%location 4
% tol = 1e-7;
% xgmax = 0.5765; % Particle Grid range
% ygmax = 0.9305; 
% xgmin = 0.576;
% ygmin = 0.930;

%location 5
% tol = 1e-7;
% xgmax = 1.314+tol; % Particle Grid range
% ygmax = 0.936+tol; 
% xgmin = 1.314-tol;
% ygmin = 0.936-tol;

%location 6
% tol = 1e-7;
% xgmax = 1.134+tol; % Particle Grid range
% ygmax = 0.918+tol; 
% xgmin = 1.13-tol;
% ygmin = 0.918-tol;

xgrid = linspace(xgmin,xgmax,3);%xgmin:dx:xgmax;
ygrid = linspace(ygmin,ygmax,3);%ygmin:dx:ygmax;
[x0,y0] = meshgrid(xgrid,ygrid);  % grid of particles for t=0

dx = abs(xgrid(2)-xgrid(1))
dy = dx;
yIC(1,:,:) = x0'; %for scatter point display
yIC(2,:,:) = y0';
if int == 'f'
    sgn = 1;
    tVec =0:(T/(dt)-1);
else
    sgn = -1;
    tVec = flip(0:(T/(dt))-1);
end

color = 'k';
yin = yIC; %Show at the start (t=0) for video purposes
yout = yin; %Use to compute other points (t>0)
intLength = (T/(dt));

xT = zeros(length(xgrid),length(ygrid),(T/(dt)));
yT = zeros(length(xgrid),length(ygrid),(T/(dt)));

%xT(:,:,1) = reshape(yout(1,:,:),length(xgrid),length(ygrid));
%yT(:,:,1) = reshape(yout(2,:,:),length(xgrid),length(ygrid));

%% Compute bound of MSE on gyre flow, sweep through variable t
n=tic;
t0=0;
r=1;
tend = intLength + t0 - 1;
mseIntegral = zeros(length(xgrid),length(ygrid),intLength);
deltaInfty = zeros(length(xgrid),length(ygrid),intLength);
for t=t0:tend 
    t 
    tic
    for s=(t0):t+1 %variable length from t0 to t-t0       
        if int == 'f'
            sgn = 1;
            sVec1 =t0:s;
            sVec2 =s:t;
        else
            sgn = -1;
            sVec1 =flip(t0:s);
            sVec2 = flip(s:t);
        end
        
        yin = yIC;
        
        %Advect from t0->s
        for i=sVec1
            time = i*dt;
            if s == t0 % special case, no need to advect
                break
            end
            yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,0,0,0,0),sgn*dt,sgn*time,yin);
            yin = yout;
        end
        %compute the amplitude value at time s
        g0(:,:,:,s+1) = gVEC(time,yin,omega,bt,at);

        %Advect from s->t
        for i=sVec2
            time = i*dt;
            yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,0,0,0,0),sgn*dt,sgn*time,yin);
            yin = yout;
        end
        %close(v)
        % reshape 3-dim array into 2-dim array
        xT = reshape(yout(1,:,:),length(xgrid),length(ygrid));
        yT = reshape(yout(2,:,:),length(xgrid),length(ygrid));
        
        %****** dx should be way smaller
        [dxTdx0,dxTdy0] = gradient(xT,dx,dx);
        [dyTdx0,dyTdy0] = gradient(yT,dx,dx);
        % Find max eigenvalue for the specific s duration
        for i=1:length(xgrid)
            for j=1:length(ygrid)
                D(1,1) = dxTdx0(i,j);
                D(1,2) = dxTdy0(i,j);
                D(2,1) = dyTdx0(i,j);
                D(2,2) = dyTdy0(i,j);
                sigma_mse(i,j,s+1) = sqrt(max(eig(D'*D))); 
                %if s == t0 % save for the zeta value (FTLE comparison)
                %   sigma_t0t(i,j) = sqrt(max(eig(D'*D))); %FTLE
                %end
            end
        end
    end %end the variable loop of s
    toc
    %for the specific end interval of integration t, compute the MSE across
    %the whole grid
    for i=1:length(xgrid)
        for j=1:length(ygrid)
            mseIntegral(i,j,r) = trapz(sigma_mse(i,j,:)*dt);
        end
    end
    gx0 = g0(1,:,:,:);
    gy0 = g0(2,:,:,:);
    deltaInfty(:,:,r) = epsn*reshape(max(sqrt(gx0.^2+gy0.^2),[],4),[length(xgrid),length(ygrid)]); %compute amplitude of mse
    clear g0 gx0 gy0 sigma_mse
    r=r+1;
    if (r > 180) && (mod(r,5) == 0)
        %save("currentboundmse_tff.mat",'deltaInfty','mseIntegral','-append')
    end
end
comptime=toc(n)

save("currentboundmse_tff_location2_b4_T24.mat")
%% Vector Field Trajectory Plus Particles Video
clear all
close all
%addpath('data')
addpath('models')
graphicsON = 0;   % flag for graphics
tstart = tic;
A = 0.1;    % parameters from Shadden 2005 Physica D
gamma = 0.25;
U=0.2;
om = 0; %2*pi/5;  % frequency of gyre oscillations
omega = 2*pi/5;
bt=4;
at = 1.4;
epsn = 0.03; %0.1; %0.02;
dt =0.025;  % timestep
dx = 0.006; %0.015; %0.005; % %
dxv = 0.0125;
T = 12*2;     % duration of integration
int = 'f'; %''f for forward integration, 'b' for backward integration;
sf = 0.4; %4;
xmax = 2; %4;  %vector range
ymax = 1; %5;
xmin = 0;
ymin = 0;
%
xvec = xmin:dxv:xmax; %vector range
yvec = ymin:dxv:ymax;
tol = 1e-10; %0.02%0.2;
%location 0
% xgmax = 0.45+tol; %1.5; %Grid range
% ygmax = 0.55+tol; %1.5;
% xgmin = 0.45-tol;
% ygmin = 0.55-tol;

%location 1
%tol = 0.10;
%xgmax = 0.48+tol; %1.5; %Grid range
%ygmax = 0.63+tol; %1.5;
%xgmin = 0.48-tol;
%ygmin = 0.63-tol;

%location 2
% tol = 0.15;
 xgmax = 1.5+tol; %1.5; %Grid range
 ygmax = 0.75+tol; %1.5;
 xgmin = 1.5-tol;
 ygmin = 0.75-tol;

%location 3
tol = 1e-7;
xgmax = 1.0265; % Particle Grid range
ygmax = 0.0665; 
xgmin = 1.026;
ygmin = 0.066;

%location 5
tol = 1e-7;
xgmax = 1.314+tol; % Particle Grid range
ygmax = 0.936+tol; 
xgmin = 1.314-tol;
ygmin = 0.936-tol;

%location 6
tol = 1e-7;
xgmax = 1.134+tol; % Particle Grid range
ygmax = 0.918+tol; 
xgmin = 1.13-tol;
ygmin = 0.918-tol;

xgrid = linspace(xgmin,xgmax,100);%xgmin:dx:xgmax;
ygrid = linspace(ygmin,ygmax,100);%ygmin:dx:ygmax;
dx = abs(xgrid(2)-xgrid(1))
dy = dx;
%[x0grid,y0grid] = meshgrid(xgrid,ygrid);  % grid of particles for later t
[x0,y0] = meshgrid(xgrid,ygrid);  % grid of particles for t=0
[x0vec,y0vec] = meshgrid(xvec,yvec);

yvIC(1,:,:) = x0vec'; %v is for vector display
yvIC(2,:,:) = y0vec';
yIC(1,:,:) = x0'; %for scatter point display
yIC(2,:,:) = y0';
if int == 'f'
    sgn = 1;
    tVec =0:(T/(dt)-1);
else
    sgn = -1;
    tVec = flip(0:(T/(dt))-1);
end

color = 'k';


%set(gcf,'Position',[100 100 600 264.15])

step=15;
color = linspace(1,10,length(xgrid)*length(ygrid));
%v = VideoWriter('test');
%v.FrameRate = 20;
%open(v)

%%
yin = yIC; %Show at the start (t=0) for video purposes
yout = yin; %Use to compute other points (t>0)
xT = zeros(length(xgrid),length(ygrid),(T/(dt)));
yT = zeros(length(xgrid),length(ygrid),(T/(dt)));

xT(:,:,1) = reshape(yout(1,:,:),length(xgrid),length(ygrid));
yT(:,:,1) = reshape(yout(2,:,:),length(xgrid),length(ygrid));
close all
if graphicsON
    set(gcf,'Position',[0 100 1000 400])
end
for i=tVec
    i
    time = i*dt;
    if graphicsON
        scatter(yout(1,:),yout(2,:),15,color,'o','filled','MarkerFaceAlpha',0.4)
        colormap(gca,"jet")
        %axis([xmin xmax ymin ymax])
        %axis([0 2 0 1])
        %xlim([0, 2])
        %ylim([0, 1])
        %set(gcf,'color','w')
        %axis equal
        %subplot(1,2,2)
        hold on
        dy = crossflowVEC(time,yvIC,A,gamma,om,epsn,omega,bt,at,U);
        quiver(yvIC(1,1:step:end,1:step:end),yvIC(2,1:step:end,1:step:end),sf*dy(1,1:step:end,1:step:end),sf*dy(2,1:step:end,1:step:end),'off','Color','#808080','linewidth',1.8);
        %axis([xmin xmax ymin ymax])
        %axis([0 2 0 1])
        %xlim([0, 2])
        %ylim([0, 1])
        set(gcf,'color','w')
        %axis equal
        box off
        drawnow
        %hold off
    end
    yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,epsn,omega,bt,at),sgn*dt,sgn*time,yin);
    yin = yout;    
    %img = getframe(1);
    %writeVideo(v,img); 
    %hold on
    %hold off
    xT(:,:,i+1) = reshape(yout(1,:,:),length(xgrid),length(ygrid));
    yT(:,:,i+1) = reshape(yout(2,:,:),length(xgrid),length(ygrid));
end
%close(v)
% reshape 3-dim array into 2-dim array
%xT = reshape(yout(1,:,:),length(xgrid),length(ygrid));
%yT = reshape(yout(2,:,:),length(xgrid),length(ygrid));

save("currentparticles_pert_location6_pert_b4_T24_point.mat")