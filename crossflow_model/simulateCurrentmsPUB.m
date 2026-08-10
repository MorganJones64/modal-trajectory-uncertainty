
%% LONG COMPUTATIONAL Compute field MSE Snapshots on gyre flow, variable t
clear all
close all
addpath('data\','models\')
A = 0.1;    % parameters from Shadden 2005 Physica D
eps = 0.03; %0.1; %0.02;
omega_gy = 0; %2*pi/5;  % frequency of gyre oscillations
omega_cross = 2*pi/5;
U = 0.2;
bt=4;
at=1.4;
gamma = 0.25;
dt =0.025;  % timestep
dx =0.006; % 0.0125; % % %
T = 15;     % duration of integration (t-t0)
intLength = T/dt; %How long the integration is for
ns = 180; %number of snapshots desired

int = 'f'; %''f for forward integration, 'b' for backward integration;
flag =0;

xgmax = 2; %4;  %vector range
ygmax = 1; %5;
xgmin = 0;
ygmin = 0;

xgrid = xgmin:dx:xgmax;
ygrid = ygmin:dx:ygmax;
[x0grid,y0grid] = meshgrid(xgrid,ygrid);  % grid of particles for later t

yIC(1,:,:) = x0grid'; %for scatter point display
yIC(2,:,:) = y0grid';

color = 'k';
yin = yIC;
r=0;

%% Compute the Model Sensitivity Field Snapshots (Haller 2020)
clc
n=tic;
phase=[1,50,100,150]
for t0=phase%0:ns
    t0
    r = t0 + 1; %for savings MSE fields
    tic
    t=intLength + t0; %end interval of integration (fixed for figure 2)
    for s=t0:t %variable length from t0 to t-t0
        if mod(s,50) == 0 %check every 50 s values the progress
            s
        end

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
            yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,omega_gy,0,0,0,0),sgn*dt,sgn*time,yin);
            yin = yout;
        end
        %compute the amplitude value at time s
        gx0(:,:,:,s+1) = gVEC(time,yin,omega_cross,bt,at);

        %Advect from s->t
        for i=sVec2
            time = i*dt;
            yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,omega_gy,0,0,0,0),sgn*dt,sgn*time,yin);
            yin = yout;
        end
        %close(v)
        % reshape 3-dim array into 2-dim array
        xT = reshape(yout(1,:,:),length(xgrid),length(ygrid));
        yT = reshape(yout(2,:,:),length(xgrid),length(ygrid));
        
        [dxTdx0,dxTdy0] = gradient(xT,dx,dx);
        [dyTdx0,dyTdy0] = gradient(yT,dx,dx);
        % Find max eigenvalue for the specific s duration
        for i=1:length(xgrid)
            for j=1:length(ygrid)
                D(1,1) = dxTdx0(i,j);
                D(1,2) = dxTdy0(i,j);
                D(2,1) = dyTdx0(i,j);
                D(2,2) = dyTdy0(i,j);
                sigma_mse(i,j,s-t0+1) = sqrt(max(eig(D'*D))); 
                if s == t0 % save for the zeta value (FTLE comparison)
                    sigma_t0t(i,j,r) = sqrt(max(eig(D'*D))); %FTLE
                end
            end
        end
    end %end the variable loop of s
    
    %for the specific end interval of integration t, compute the MSE across
    %the whole grid
    for i=1:length(xgrid)
        for j=1:length(ygrid)
            mseIntegral(i,j,r) = trapz(sigma_mse(i,j,:) * dt);
        end
    end
    deltaInfty(:,:,r) = eps*reshape(max(max(abs(gx0),[],4),[],1),[length(xgrid),length(ygrid)]); %compute amplitude of mse
    MSE(:,:,r) = (mseIntegral(:,:,r).*deltaInfty(:,:,r)).^2;
    zeta(:,:,r) = (1/T)*(log(mseIntegral(:,:,r)./sigma_t0t(:,:,r)) + log(deltaInfty(:,:,r))); %Part of equation (V.I6) (Kaszas and Haller 2020)
    toc
    if mod(r,20) == 0 %check every 50 s values the progress
        break
    end
end
comptime=toc(n)

save("data/currentMTUsnapshots.mat")
