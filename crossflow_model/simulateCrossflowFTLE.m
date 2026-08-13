%% ========================================================================
%  Compute Forward FTLE for Crossflow Model
%  ========================================================================
%  This script computes the forward FTLE field at phase 100 (half period)
%  for the kinematic crossflow model with oscillating jet perturbation.
%  
%  System: Base gyre flow + oscillating crossflow jet
%  Phase 100 corresponds to T_p/2 where T_p = 5 seconds (half oscillation)
%  ========================================================================

clear all
close all

% Add model functions
addpath('models');
addpath('data');

fprintf('\n========================================\n');
fprintf('Forward FTLE Computation - Crossflow Model\n');
fprintf('========================================\n\n');

%% ========================================================================
%  Model Parameters
%  ========================================================================

% Base gyre flow parameters (from Shadden 2005)
A = 0.1;         % Gyre amplitude
gamma = 0.25;    % Asymmetry parameter
U = 0.2;         % Background flow velocity
om = 0;          % Gyre oscillation frequency (0 = steady gyres)

% Crossflow jet parameters
omega = 2*pi/5;  % Jet oscillation frequency (period T_p = 5 seconds)
bt = 4;          % Jet width parameter
at = 1.4;        % Jet center position
eps = 0.03;      % Jet amplitude

% Numerical parameters
dt = 0.025;      % Time step
dx = 0.006;      % Spatial resolution for FTLE grid
T = 15;          % Integration duration

% Integration direction
int = 'f';       % 'f' for forward, 'b' for backward

% Domain bounds
xmin = 0;
xmax = 2;
ymin = 0;
ymax = 1;

% FTLE grid
xgrid = xmin:dx:xmax;
ygrid = ymin:dx:ymax;
[x0, y0] = meshgrid(xgrid, ygrid);

% Initial condition array
yIC(1, :, :) = x0';
yIC(2, :, :) = y0';

% Phase parameter (half period: T_p/2 = 5/2 = 2.5 sec = phase 100)
phase = 100;  % Phase 100 = t0 = 100 * dt = 2.5 seconds = T_p/2

fprintf('Model Parameters:\n');
fprintf('  Base flow amplitude (A): %.2f\n', A);
fprintf('  Background velocity (U): %.2f\n', U);
fprintf('  Jet amplitude (eps): %.3f\n', eps);
fprintf('  Jet frequency (omega): %.4f rad/s (T_p = %.2f s)\n', omega, 2*pi/omega);
fprintf('  Integration time (T): %.1f seconds\n', T);
fprintf('  Time step (dt): %.4f seconds\n', dt);
fprintf('  Grid resolution (dx): %.4f\n', dx);
fprintf('  Phase: %d (t0 = %.2f s = T_p/2)\n\n', phase, phase*dt);


%% ========================================================================
%  Integrate Particle Trajectories
%  ========================================================================

fprintf('Integrating particle trajectories...\n');

% Set integration direction and time vector
if strcmp(int, 'f')
    sgn = 1;
    tVec = phase:(T/dt + phase);
else
    sgn = -1;
    tVec = flip(phase:(T/dt + phase));
end

yin = yIC;

% Integrate trajectories
tic;
for i = tVec
    time = i * dt;
    yout = rk4singlestep(@(t,x) crossflowVEC(t, x, U, A, gamma, om, eps, omega, bt, at), ...
        sgn*dt, sgn*time, yin);
    yin = yout;
end
elapsed = toc;

fprintf('Integration complete! Time: %.2f seconds\n\n', elapsed);

% Reshape to 2D arrays
xT = reshape(yout(1, :, :), length(xgrid), length(ygrid));
yT = reshape(yout(2, :, :), length(xgrid), length(ygrid));


%% ========================================================================
%  Compute FTLE Field
%  ========================================================================

fprintf('Computing FTLE field...\n');

% Compute gradients using finite differences
[dxTdx0, dxTdy0] = gradient(xT, dx, dx);
[dyTdx0, dyTdy0] = gradient(yT, dx, dx);

% Compute FTLE at each grid point
for i = 1:length(xgrid)
    for j = 1:length(ygrid)
        % Deformation gradient tensor
        D = [dxTdx0(i,j), dxTdy0(i,j); ...
             dyTdx0(i,j), dyTdy0(i,j)];
        
        % FTLE: (1/T) * log(sqrt(lambda_max))
        sigma(i, j) = (1/T) * log(sqrt(max(eig(D'*D))));
    end
end

fprintf('FTLE computation complete!\n\n');

% Save results
filename = sprintf('FTLE_crossflow_phase%d.mat', phase);
save(filename, 'sigma', 'x0', 'y0', 'xgrid', 'ygrid', 'T', 'dt', 'phase', ...
    'A', 'U', 'eps', 'omega', 'bt', 'at');
fprintf('Results saved to: %s\n\n', filename);

%% ========================================================================
%  Visualize FTLE Field
%  ========================================================================

fprintf('Generating visualization...\n');

% Set colormap based on integration direction
if strcmp(int, 'f')
    mycolormap = customcolormap(linspace(0,1,7),{'#FADA0E','#EE2020','#000000','#000000','#000000','#000000','#000000'}); %red
    titleStr = 'Forward FTLE';
else
    mycolormap = customcolormap(linspace(0,1,6),{'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'});
    titleStr = 'Backward FTLE';
end

% Create figure
figure('Position', [100, 100, 1000, 500]);
contourf(x0', y0', sigma, 200, 'LineStyle', 'none');
axis([xmin, xmax, ymin, ymax]);
axis equal tight;

% Formatting
colormap(mycolormap);
clim([0, 0.25]);
c = colorbar;
c.FontSize = 20;
c.Ticks = [0, 0.125, 0.25];
c.FontName = 'Times New Roman';

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
title(sprintf('%s (Phase %d, $t_0 = %.2f$ s)', titleStr, phase, phase*dt), ...
    'Interpreter', 'latex', 'FontSize', 22);

set(gcf, 'color', 'w');
box on;

fprintf('\n========================================\n');
fprintf('FTLE Computation Complete!\n');
fprintf('========================================\n\n');


