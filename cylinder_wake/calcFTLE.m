%% ========================================================================
%  Compute Backward FTLE Fields for Cylinder Wake
%  ========================================================================
%  This script computes backward FTLE fields for:
%  1. Original flow field (u, v from cylindervelocity_zen.mat)
%  2. Modal representation (umean + u1 + u2 from cylinderPODmodesPair.mat)
%  
%  Based on the paper (main.tex), equation (cyl-sys2):
%  tilde{u} = umean + u1 + u2 (modal representation)
%  u' = u3 (perturbation)
%  ========================================================================

clear all

% Add necessary paths
addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

% Load cylinder coordinates for plotting
load('cylcoords.mat');

%% ========================================================================
%  Part 1: Compute FTLE for Original Flow Field
%  ========================================================================

fprintf('\n========================================\n');
fprintf('PART 1: Computing FTLE for Original Flow\n');
fprintf('========================================\n\n');

% Load original velocity data
load('cylindervelocity_zen.mat');

% Set velocity fields
uMesh = u;
vMesh = v;

% FTLE computation parameters
xMinROI = -2;
xMaxROI = 20;
yMinROI = -5;
yMaxROI = 5;
ROInx = 800;
ROIny = 400;
method = 'RK4';

% Time parameters
tVec = (0:1:(nt-1)) * dt';
yVec = y;
xVec = x;

% Backward FTLE parameters
tLength = -199;  % Negative for backward integration
tStep = -1;      % Negative for backward integration
fstart = 200;    % Starting frame

fprintf('FTLE Parameters:\n');
fprintf('  Integration time: %d frames\n', abs(tLength));
fprintf('  Starting frame: %d\n', fstart);
fprintf('  ROI: x=[%.1f, %.1f], y=[%.1f, %.1f]\n', xMinROI, xMaxROI, yMinROI, yMaxROI);
fprintf('  Grid resolution: %d x %d\n', ROInx, ROIny);
fprintf('  Integration method: %s\n\n', method);

% Create extrapolation functions for FTLE
fprintf('Creating interpolation functions...\n');
Fx = interpolantxy(uMesh, xVec, yVec, tVec, 1);
Fy = interpolantxy(vMesh, xVec, yVec, tVec, 0);

uExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx);
vExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy);

% Compute FTLE
fprintf('Computing backward FTLE for original flow...\n');
[sigma, xPos, yPos] = FTLE(uMesh, vMesh, ...
    xVec, yVec, ...
    fstart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ...
    ROInx, ROIny, method, ...
    'extrap', true, 'uExtrap', uExtrap, 'vExtrap', vExtrap);

% Save results
fprintf('Saving results to bFTLE_original.mat...\n');
save('bFTLE_original.mat', 'sigma', 'xPos', 'yPos', ...
    'tLength', 'tVec', 'dt', 'xMaxROI', 'xMinROI', 'yMaxROI', 'yMinROI', ...
    'fstart', 'method');

fprintf('Original flow FTLE computation complete!\n\n');

%% ========================================================================
%  Part 2: Compute FTLE for Modal Representation (umean + u1 + u2)
%  ========================================================================

fprintf('========================================\n');
fprintf('PART 2: Computing FTLE for Modal Representation\n');
fprintf('========================================\n\n');

% Load POD modes
load('cylinderPODmodesPair.mat');

% Construct modal representation: tilde{u} = umean + u1 + u2
fprintf('Constructing modal representation: tilde{u} = umean + u1 + u2\n');
uMesh_modal = umean + u1 + u2;
vMesh_modal = vmean + v1 + v2;

% Create extrapolation functions for modal representation
fprintf('Creating interpolation functions for modal representation...\n');
Fx_modal = interpolantxy(uMesh_modal, xVec, yVec, tVec, 1);
Fy_modal = interpolantxy(vMesh_modal, xVec, yVec, tVec, 0);

uExtrap_modal = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx_modal);
vExtrap_modal = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy_modal);

% Compute FTLE
fprintf('Computing backward FTLE for modal representation...\n');
[sigma_modal, xPos_modal, yPos_modal] = FTLE(uMesh_modal, vMesh_modal, ...
    xVec, yVec, ...
    fstart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ...
    ROInx, ROIny, method, ...
    'extrap', true, 'uExtrap', uExtrap_modal, 'vExtrap', vExtrap_modal);

% Save results
fprintf('Saving results to bFTLE_modalrep.mat...\n');
save('bFTLE_modalrep.mat', 'sigma_modal', 'xPos_modal', 'yPos_modal', ...
    'tLength', 'tVec', 'dt', 'xMaxROI', 'xMinROI', 'yMaxROI', 'yMinROI', ...
    'fstart', 'method');

fprintf('Modal representation FTLE computation complete!\n\n');

%% ========================================================================
%  Part 3: Visualize FTLE Fields
%  ========================================================================

fprintf('========================================\n');
fprintf('PART 3: Visualizing FTLE Fields\n');
fprintf('========================================\n\n');

% Compute FTLE from sigma values
FTLE_original = (1 / abs(tLength * dt)) * log(sigma);
FTLE_modal = (1 / abs(tLength * dt)) * log(sigma_modal);

% Get grid coordinates
xg = xPos(:, :, 1);
yg = yPos(:, :, 1);

% Plotting parameters
cmap = customcolormap(linspace(0, 1, 5), {'#54ff82', '#1fe9ff', '#1FB8FF', '#184ADC', '#000000'});
width = 800;
height = 320;
axislength = 20;
ep = 0.01;
umin = 0;
umax = 0.25;
SetFont = 'Times New Roman';

% Plot Original FTLE
fprintf('Plotting original FTLE field...\n');
figure(1);
set(gcf, 'Position', [100, 100, width, height]);
contourf(xg, yg, FTLE_original, 200, 'LineStyle', 'none');
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
title('Backward FTLE: Original Flow', 'Interpreter', 'latex', 'fontsize', axislength);
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

% Plot Modal Representation FTLE
fprintf('Plotting modal representation FTLE field...\n');
figure(2);
set(gcf, 'Position', [100, 500, width, height]);
contourf(xg, yg, FTLE_modal, 200, 'LineStyle', 'none');
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
title('Backward FTLE: Modal Representation ($\tilde{u} = \bar{u} + u_1 + u_2$)', ...
    'Interpreter', 'latex', 'fontsize', axislength);
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
fprintf('All FTLE computations and visualizations complete!\n');
fprintf('========================================\n');
fprintf('\nOutput files:\n');
fprintf('  - bFTLE_original.mat (original flow field)\n');
fprintf('  - bFTLE_modalrep.mat (modal representation: umean + u1 + u2)\n');
fprintf('\nFigures:\n');
fprintf('  - Figure 1: Original flow FTLE\n');
fprintf('  - Figure 2: Modal representation FTLE\n\n');
%%
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

