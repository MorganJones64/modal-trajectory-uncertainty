
%% ========================================================================
%  Compute Modal-Trajectory Uncertainty for Crossflow Model
%  ========================================================================
%  This script computes MTU at phase 100 (half period) for the kinematic
%  crossflow model. Compares base gyre flow vs. gyre + oscillating jet.
%  
%  Modal representation: Base gyre flow (no jet)
%  Perturbation: Oscillating crossflow jet
%  ========================================================================

clear all
close all
addpath('data/','models/')

fprintf('\n========================================\n');
fprintf('MTU Computation - Crossflow Model\n');
fprintf('========================================\n\n');

%========================== Model Parameters ===========================
A = 0.1;            % Gyre amplitude
gamma = 0.25;       % Asymmetry parameter 
U = 0.2;            % Background velocity
omega_gy = 0;       % Gyre oscillation frequency (steady)

eps = 0.03;         % Jet amplitude
omega_cross = 2*pi/5; % Jet frequency (T_p = 5 s)
bt = 4;             % Jet width
at = 1.4;           % Jet center

dt = 0.025;         % Time step
dx = 0.006;         % Spatial resolution
T = 15;             % Integration duration
intLength = T/dt;

int = 'f';          % Forward integration
phase = 100;        % t0 = 2.5 s = T_p/2

% Domain
xgmin = 0;
xgmax = 2;
ygmin = 0;
ygmax = 1;

fprintf('Model Parameters:\n');
fprintf('  Base flow: A=%.2f, U=%.2f\n', A, U);
fprintf('  Jet: eps=%.3f, omega=%.4f\n', eps, omega_cross);
fprintf('  Phase: %d (t0=%.2f s)\n\n', phase, phase*dt);

% Setup grid
xgrid = xgmin:dx:xgmax;
ygrid = ygmin:dx:ygmax;
[x0grid,y0grid] = meshgrid(xgrid,ygrid);

yIC(1,:,:) = x0grid';
yIC(2,:,:) = y0grid';

%% Compute MTU
fprintf('Computing MTU at phase %d...\n', phase);

t0 = phase;
t = intLength + t0;

if strcmp(int, 'f')
    sgn = 1;
else
    sgn = -1;
end

num_times = t - t0 + 1;
sigma_mse = zeros(length(xgrid), length(ygrid), num_times);
g0 = zeros(2, length(xgrid), length(ygrid), num_times);

tic;
for s=t0:t
            if mod(s-t0,50)==0
        fprintf('  Progress: %.1f%%\n', 100*(s-t0)/(t-t0));
    end
    
    if strcmp(int,'f')
        sVec1 = t0:s;
        sVec2 = s:t;
    else
        sVec1 = flip(t0:s);
        sVec2 = flip(s:t);
    end
    
    yin = yIC;
    
    % Advect t0->s (base flow only)
    if s > t0
        for i=sVec1
            time = i*dt;
            yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,omega_gy,0,0,0,0),sgn*dt,sgn*time,yin);
            yin = yout;
        end
    else
        time = t0*dt; %special case: s=t0, no advection needed
    end
    
    % Save perturbation amplitude at s
    g0(:,:,:,s+1) = gVEC(time,yin,omega_cross,bt,at);
    
    % Advect s->t (base flow only)
    for i=sVec2
        time = i*dt;
        yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,omega_gy,0,0,0,0),sgn*dt,sgn*time,yin);
        yin = yout;
    end
    
    % Reshape
    xT = reshape(yout(1,:,:),length(xgrid),length(ygrid));
    yT = reshape(yout(2,:,:),length(xgrid),length(ygrid));
    
    % Gradients
    [dxTdx0,dxTdy0] = gradient(xT,dx,dx);
    [dyTdx0,dyTdy0] = gradient(yT,dx,dx);
    
    % Compute eigenvalues
    for i=1:length(xgrid)
        for j=1:length(ygrid)
            D = [dxTdx0(i,j), dxTdy0(i,j); dyTdx0(i,j), dyTdy0(i,j)];
            sigma_mse(i,j,s-t0+1) = sqrt(max(eig(D'*D)));
            if s==t0
                sigma_ftle(i,j) = sqrt(max(eig(D'*D)));
            end
        end
    end
end
comptime=toc;

fprintf('\nComputation complete! Time: %.2f s\n\n', comptime);

%======================= Compute MTU Metrics ===========================
fprintf('Computing MTU metrics...\n');

mseIntegral = zeros(length(xgrid),length(ygrid));
for i=1:length(xgrid)
    for j=1:length(ygrid)
        mseIntegral(i,j) = trapz(squeeze(sigma_mse(i,j,:))*dt);
    end
end

gx0 = g0(1,:,:,:);
gy0 = g0(2,:,:,:);
deltaInfty = eps*reshape(max(sqrt(gx0.^2+gy0.^2),[],4), [length(xgrid),length(ygrid)]);
scaledMTU = log((mseIntegral.*deltaInfty).^2) / (1 / (2*T));
zeta = (1/T)*(log(mseIntegral./sigma_ftle) + log(deltaInfty));

fprintf('MTU metrics computed!\n\n');

% Save Results
filename = sprintf('data/MTU_crossflow_phase%d.mat', phase);
save(filename, 'scaledMTU', 'zeta', 'mseIntegral', 'deltaInfty', 'sigma_ftle', ...
    'x0grid', 'y0grid', 'xgrid', 'ygrid', 'phase', 'T', 'dt', 'dx', ...
    'A', 'U', 'gamma', 'eps', 'omega_cross', 'bt', 'at', 'comptime');
fprintf('Saved: %s\n\n', filename);

%======================= Visualize =====================================
fprintf('Generating figure...\n');
umin = -0.1;
umax = 0.15;
mycolormap = customcolormap(linspace(0,1,7),{'#FADA0E','#EE2020','#000000','#000000','#000000','#000000','#000000'}); %red
figure('Position',[100,100,900,400]);
contourf(x0grid,y0grid,scaledMTU(:,:,1)',200,'LineStyle','none');
axis equal tight; xlim([xgmin,xgmax]); ylim([ygmin,ygmax]);
colormap(mycolormap); colorbar;
xlabel('$x$','Interpreter','latex','FontSize',18);
ylabel('$y$','Interpreter','latex','FontSize',18);
clim([umin, umax]);
title(sprintf('Scaled MTU (Phase %d)',phase),'Interpreter','latex','FontSize',20);
set(gcf,'color','w'); box on;
