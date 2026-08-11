%% ========================================================================
%  Compute Modal-Trajectory Uncertainty (MTU) for Cylinder Wake
%  ========================================================================
%  This script computes backward modal-trajectory uncertainty (MTU) for:
%  Modal representation: tilde{u} = umean + u1 + u2
%  Perturbation: u' = u3
%  
%  Based on the paper (main.tex), this quantifies how mode u3 influences
%  the Lagrangian coherent structures of the baseline system tilde{u}.
%  
%  This version is streamlined for backward FTLE computation at a single frame.
%  For a general case (forward/backward, multiple frames), see calcMTU_general.m
%  ========================================================================

clear all

% Add necessary paths
addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

% Load cylinder coordinates for plotting
load('cylcoords.mat');

fprintf('\n========================================\n');
fprintf('Modal-Trajectory Uncertainty Computation\n');
fprintf('========================================\n\n');

%% ========================================================================
%  Load Data and Construct Modal Fields
%  ========================================================================

fprintf('Loading data...\n');
load('cylindervelocity_zen.mat');
load('cylinderPODmodesPair.mat');

% Construct modal representation and perturbation
fprintf('\nConstructing modal fields:\n');
fprintf('  Modal representation: tilde{u} = umean + u1 + u2\n');
fprintf('  Perturbation: u'' = u3\n\n');

ufMesh = umean + u1 + u2;  % Modal representation (baseline)
vfMesh = vmean + v1 + v2;
ugMesh = u3;               % Perturbation
vgMesh = v3;

% Clear unnecessary variables
clear u1 u2 u3 u4 ut umean v1 v2 v3 v4 vt vmean

%% ========================================================================
%  Setup MTU Computation Parameters
%  ========================================================================

% Region of interest
xMinROI = -2;
xMaxROI = 20;
yMinROI = -5;
yMaxROI = 5;
ROInx = 800;
ROIny = 400;

% Integration method
method = 'RK4';

% Time parameters
tVec = (0:1:(nt-1)) * dt';
yVec = y;
xVec = x;

% Backward FTLE parameters (single frame)
tLength = -199;  % Negative for backward integration
tStep = -1;      % Negative for backward integration
fstart = 200;    % Starting frame

fprintf('MTU Parameters:\n');
fprintf('  Integration time: %d frames\n', abs(tLength));
fprintf('  Starting frame: %d\n', fstart);
fprintf('  ROI: x=[%.1f, %.1f], y=[%.1f, %.1f]\n', xMinROI, xMaxROI, yMinROI, yMaxROI);
fprintf('  Grid resolution: %d x %d\n', ROInx, ROIny);
fprintf('  Integration method: %s\n\n', method);

%% ========================================================================
%  Create Interpolation Functions
%  ========================================================================

fprintf('Creating interpolation functions...\n');
Ffx = interpolantxy(ufMesh, xVec, yVec, tVec, 1);
Ffy = interpolantxy(vfMesh, xVec, yVec, tVec, 0);
Fgx = interpolantxy(ugMesh, xVec, yVec, tVec, 0);
Fgy = interpolantxy(vgMesh, xVec, yVec, tVec, 0);

ufExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Ffx);
vfExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Ffy);
ugExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fgx);
vgExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fgy);
clear Ffx Ffy Fgx Fgy

%% ========================================================================
%  Compute MTU
%  ========================================================================

fprintf('\nComputing modal-trajectory uncertainty...\n');
tic;
[sigma_ftle, cseIntegral, deltaInfty, xPos, yPos] = MTU(ufMesh, vfMesh, ...
    ugMesh, vgMesh, ...
    xVec, yVec, ...
    fstart, tLength, tStep, fstart, eps, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ...
    ROInx, ROIny, method, ...
    'extrap', true, 'ufExtrap', ufExtrap, 'vfExtrap', vfExtrap, ...
    'ugExtrap', ugExtrap, 'vgExtrap', vgExtrap);
toc;

% Save results
filename = 'MTU_cylinder_wake.mat';
notes = 'Modal representation: f=umean+u1+u2, Perturbation: g=u3, backward MTU';
fprintf('\nSaving results to %s...\n', filename);
save(filename, 'sigma_ftle', 'cseIntegral', 'deltaInfty', 'xPos', 'yPos', ...
    'dt', 'tLength', 'xMaxROI', 'xMinROI', 'yMaxROI', 'yMinROI', ...
    'fstart', 'method', 'notes');

fprintf('MTU computation complete!\n\n');

%% ========================================================================
%  Compute MTU Metrics
%  ========================================================================

fprintf('\n========================================\n');
fprintf('Computing MTU Metrics\n');
fprintf('========================================\n\n');

% Compute FTLE
FTLE = (1 / abs(tLength * dt)) * log(sigma_ftle);

% Compute scaled MTU (modal-trajectory uncertainty)
scaledMTU = log((deltaInfty .* cseIntegral).^2) / abs(2 * tLength * dt);

% Compute FTLE perturbation (zeta)
zeta = (1 / abs(tLength * dt)) * (log(deltaInfty) + log(cseIntegral) - log(sigma_ftle));

% Get grid coordinates
x = xPos(:, :, 1);
y = yPos(:, :, 1);

fprintf('Metrics computed successfully!\n\n');
%% ========================================================================
%  Visualize Scaled MTU Field
%  ========================================================================

fprintf('========================================\n');
fprintf('Visualizing Scaled MTU Field\n');
fprintf('========================================\n\n');

% Plotting parameters
cmap = customcolormap(linspace(0, 1, 6), {'#54ff82', '#1fe9ff', '#1FB8FF', '#184ADC', '#000000', '#000000'});
width = 800;
height = 320;
axislength = 20;
ep = 0.01;
SetFont = 'Times New Roman';

% Color limits for scaled MTU
umin = -0.45;
umax = -0.25;

% Plot scaled MTU field
figure(1);
set(gcf, 'Position', [100, 100, width, height]);
contourf(x, y, scaledMTU, 200, 'LineStyle', 'none');
hold on;
colormap(cmap);
axis equal;
clim([umin, umax]);

c = colorbar;
c.Ticks = linspace(umin, umax, 3);
c.FontSize = axislength;
c.FontName = SetFont;
c.LineWidth = 1.2;
c.TickLength = 0.02;

xlabel('$x$', 'Interpreter', 'latex', 'fontsize', axislength);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', axislength);
title('Scaled Modal-Trajectory Uncertainty', 'Interpreter', 'latex', 'fontsize', axislength);
set(gca, 'LineWidth', 1.2, 'TickLength', [ep, ep]);
set(gcf, 'color', 'w');

% Add cylinder
fill(xcyl, ycyl, [1, 1, 1]);
plot(xcyl, ycyl, 'k', 'LineWidth', 1.2);

xlim([-2, 20]);
ylim([-4, 4]);
xticks(0:4:20);
yticks(-4:4:4);

ax = gca;
box on;
ax.FontName = SetFont;
ax.FontSize = axislength;
hold off;
drawnow;

fprintf('\n========================================\n');
fprintf('MTU Analysis Complete!\n');
fprintf('========================================\n');
fprintf('\nOutput files:\n');
fprintf('  - %s\n', filename);
fprintf('\nFigures:\n');
fprintf('  - Figure 1: Scaled MTU field\n\n');

%% ========================================================================
%  Helper Functions
%  ========================================================================

function F = interpolantxy(uvec,xVec,yVec,tVec,bcs)
    nt = length(tVec);
    nx = length(xVec);
    ny = length(yVec);
    yb = 1;
    xb = 1;
    xVecExt = [xVec(1)-xb, xVec];
    yVecExt = [yVec(1)-yb, yVec, yVec(end)+yb];
    if bcs == 0
        uleft = zeros(ny,1);
        utb = zeros(1,nx+1);
    else
        uleft = ones(ny,1);
        utb = ones(1,nx+1);
    end
    for i = 1:nt
        utemp = [uleft, uvec(:,:,i)]; %fix
        uExt(:,:,i) = [utb; utemp; utb];
    end
    [yi,xi,ti] = ndgrid(yVecExt,xVecExt, tVec);
    F = griddedInterpolant(yi, xi, ti, uExt,'makima','nearest'); 
end

function u = uExtrapolateM(x, y, t, u, F) 
    u = F(y,x,t);
end

