%% ========================================================================
%  Complete Point-Wise MTU Analysis for Cylinder Wake (Figure 7)
%  ========================================================================
%  This script performs the complete analysis for Figure 7 in the paper:
%  1. Computes MTU bounds at probe points
%  2. Computes particle trajectories for original flow
%  3. Computes particle trajectories for modal representation
%  4. Plots upper bounds vs actual trajectory errors
%  
%  Modal representation: tilde{u} = umean + u1 + u2
%  Perturbation: u' = u3
%  ========================================================================

clear all

% Add necessary paths
addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

fprintf('\n========================================\n');
fprintf('Point-Wise MTU Bound Computation\n');
fprintf('========================================\n\n');

%% ========================================================================
%  Load Data and Construct Modal Fields
%  ========================================================================

fprintf('Loading data...\n');
load('cylindervelocity_zen.mat');
load('cylinderPODmodesPair.mat');

% Construct modal representation and perturbation
ufMesh = umean + u1 + u2;  % Modal representation
vfMesh = vmean + v1 + v2;
ugMesh = u3;               % Perturbation
vgMesh = v3;

% Clear unnecessary variables
clear u1 u2 u3 u4 ut umean v1 v2 v3 v4 vt vmean

%% ========================================================================
%  Setup Computation Parameters
%  ========================================================================

% Time parameters
tVec = (0:1:(nt-1)) * dt';
yVec = y;
xVec = x;

% Create interpolation functions
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
%  Define Probe Points
%  ========================================================================

% Probe point locations (Figure 7 in paper)
% Points are selected along the vortex street at different streamwise positions
pxList = [15.68, 11.189, 4.33];  % Streamwise positions
pyList = [1.32, 1.942, 0.04];    % Lateral positions
pNames = {'p11', 'p13', 'p6c'};  % Point names for file saving

fprintf('\nProbe points:\n');
for i = 1:length(pxList)
    fprintf('  %s: (x, y) = (%.3f, %.3f)\n', pNames{i}, pxList(i), pyList(i));
end
fprintf('\n');

% ROI parameters for each point
dtol = 1e-3;  % Small tolerance around each point
ROInx = 3;    % 3x3 grid around each point
ROIny = 3;
method = 'RK4';

% Backward FTLE parameters
tLength = -199;  % Negative for backward integration
finc = 1;
tStep = -1;
fstart = 200;

%% ========================================================================
%  Loop Over All Probe Points
%  ========================================================================

for p = 1:length(pxList)
    
    % Current point
    px = pxList(p);
    py = pyList(p);
    
    fprintf('\n========================================\n');
    fprintf('Processing point %d/%d: %s\n', p, length(pxList), pNames{p});
    fprintf('Location: (x, y) = (%.3f, %.3f)\n', px, py);
    fprintf('========================================\n\n');
    
    % Define ROI around point
    xMinROI = px - dtol;
    xMaxROI = px + dtol;
    yMinROI = py - dtol;
    yMaxROI = py + dtol;
    
    % Preallocate arrays
    sigma_ftle = zeros([ROInx, ROIny, abs(tLength)]);
    cseIntegral = zeros([ROInx, ROIny, abs(tLength)]);
    deltaInfty = zeros([ROInx, ROIny, abs(tLength)]);
    
    % Metadata
    notes = sprintf('Baseline: umean+u1+u2, Perturbation: u3, backward bound, %s', pNames{p});
    filename = sprintf('bounddatasave_backwards_%s_fmnm1m2_gm3.mat', pNames{p});
    
    % Initial save with parameters
    save(filename, ...
        'sigma_ftle', 'cseIntegral', 'deltaInfty', ...
        'dt', 'tLength', ...
        'xMaxROI', 'xMinROI', 'yMaxROI', 'yMinROI', ...
        'fstart', 'method', 'notes');
    
    % Run MTU bound computation
    fprintf('Computing MTU bounds...\n');
    tic;
    
    [cseIntegral, deltaInfty, xPos, yPos] = MTUbound( ...
        ufMesh, vfMesh, ...
        ugMesh, vgMesh, ...
        xVec, yVec, ...
        fstart, tLength, tStep, dt, filename, ...
        xMinROI, xMaxROI, yMinROI, yMaxROI, ...
        ROInx, ROIny, method, ...
        'extrap', true, ...
        'ufExtrap', ufExtrap, ...
        'vfExtrap', vfExtrap, ...
        'ugExtrap', ugExtrap, ...
        'vgExtrap', vgExtrap);
    
    elapsed = toc;
    fprintf('Computation time: %.2f seconds\n', elapsed);
    
    % Append results
    save(filename, 'cseIntegral', 'deltaInfty', 'xPos', 'yPos', '-append');
    fprintf('Saved: %s\n', filename);
    
end

fprintf('\n========================================\n');
fprintf('All point-wise MTU computations complete!\n');
fprintf('========================================\n\n');

%% ========================================================================
%  Part 2: Compute Particle Positions for Original Flow
%  ========================================================================

fprintf('\n========================================\n');
fprintf('PART 2: Particle Position Computation\n');
fprintf('========================================\n\n');

% Reload original flow data
fprintf('Loading original flow data...\n');
load('cylindervelocity_zen.mat');

% Use original flow field
uMesh_orig = u;
vMesh_orig = v;

% Create extrapolation functions for original flow
fprintf('Creating interpolation functions for original flow...\n');
Fx_orig = interpolantxy(uMesh_orig, xVec, yVec, tVec, 1);
Fy_orig = interpolantxy(vMesh_orig, xVec, yVec, tVec, 0);

uExtrap_orig = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx_orig);
vExtrap_orig = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy_orig);

% Loop over probe points for original flow
for p = 1:length(pxList)
    
    px = pxList(p);
    py = pyList(p);
    
    fprintf('\n========================================\n');
    fprintf('Processing point %d/%d: %s (Original Flow)\n', p, length(pxList), pNames{p});
    fprintf('Location: (x, y) = (%.3f, %.3f)\n', px, py);
    fprintf('========================================\n\n');
    
    % Define ROI
    xMinROI = px - dtol;
    xMaxROI = px + dtol;
    yMinROI = py - dtol;
    yMaxROI = py + dtol;
    
    % Output filename
    filename_orig = sprintf('particlePos_backwards_%s_original.mat', pNames{p});
    
    % Compute particle trajectories
    fprintf('Computing particle trajectories for original flow...\n');
    tic;
    
    [sigma_ftle_orig, xPos_orig, yPos_orig] = FTLEsb( ...
        uMesh_orig, vMesh_orig, ...
        xVec, yVec, ...
        fstart, tLength, tStep, dt, ...
        xMinROI, xMaxROI, ...
        yMinROI, yMaxROI, ...
        ROInx, ROIny, method, ...
        'extrap', true, ...
        'uExtrap', uExtrap_orig, ...
        'vExtrap', vExtrap_orig);
    
    elapsed = toc;
    fprintf('Computation time: %.2f seconds\n', elapsed);
    
    % Save results
    save(filename_orig, ...
        'xPos', 'yPos', 'sigma_ftle', ...
        'yMaxROI', 'yMinROI', 'xMaxROI', 'xMinROI', ...
        'tLength', 'tVec', 'dt', 'fstart', 'finc', 'dtol');
    
    % Store variables with original suffix for later use
    eval(sprintf('xPos_orig_%s = xPos_orig;', pNames{p}));
    eval(sprintf('yPos_orig_%s = yPos_orig;', pNames{p}));
    
    fprintf('Saved: %s\n', filename_orig);
    
end

clear uMesh_orig vMesh_orig Fx_orig Fy_orig uExtrap_orig vExtrap_orig

%% ========================================================================
%  Part 3: Compute Particle Positions for Modal Representation
%  ========================================================================

fprintf('\n========================================\n');
fprintf('PART 3: Particle Positions (Modal Rep)\n');
fprintf('========================================\n\n');

% Reload POD modes
fprintf('Loading POD modes...\n');
load('cylinderPODmodesPair.mat');

% Construct modal representation (same as in Part 1)
uMesh_modal = umean + u1 + u2;
vMesh_modal = vmean + v1 + v2;

% Create extrapolation functions for modal representation
fprintf('Creating interpolation functions for modal representation...\n');
Fx_modal = interpolantxy(uMesh_modal, xVec, yVec, tVec, 1);
Fy_modal = interpolantxy(vMesh_modal, xVec, yVec, tVec, 0);

uExtrap_modal = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx_modal);
vExtrap_modal = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy_modal);

% Loop over probe points for modal representation
for p = 1:length(pxList)
    
    px = pxList(p);
    py = pyList(p);
    
    fprintf('\n========================================\n');
    fprintf('Processing point %d/%d: %s (Modal Rep)\n', p, length(pxList), pNames{p});
    fprintf('Location: (x, y) = (%.3f, %.3f)\n', px, py);
    fprintf('========================================\n\n');
    
    % Define ROI
    xMinROI = px - dtol;
    xMaxROI = px + dtol;
    yMinROI = py - dtol;
    yMaxROI = py + dtol;
    
    % Output filename
    filename_modal = sprintf('particlePos_backwards_%s_fmnm1m2.mat', pNames{p});
    
    % Compute particle trajectories
    fprintf('Computing particle trajectories for modal representation...\n');
    tic;
    
    [sigma_ftle_modal, xPos_modal, yPos_modal] = FTLEsb( ...
        uMesh_modal, vMesh_modal, ...
        xVec, yVec, ...
        fstart, tLength, tStep, dt, ...
        xMinROI, xMaxROI, ...
        yMinROI, yMaxROI, ...
        ROInx, ROIny, method, ...
        'extrap', true, ...
        'uExtrap', uExtrap_modal, ...
        'vExtrap', vExtrap_modal);
    
    elapsed = toc;
    fprintf('Computation time: %.2f seconds\n', elapsed);
    
    % Save results
    save(filename_modal, ...
        'xPos', 'yPos', 'sigma_ftle', ...
        'yMaxROI', 'yMinROI', 'xMaxROI', 'xMinROI', ...
        'tLength', 'tVec', 'dt', 'fstart', 'finc', 'dtol');
    
    % Store variables with modal suffix for later use
    eval(sprintf('xPos_modal_%s = xPos_modal;', pNames{p}));
    eval(sprintf('yPos_modal_%s = yPos_modal;', pNames{p}));
    
    fprintf('Saved: %s\n', filename_modal);
    
end

clear uMesh_modal vMesh_modal Fx_modal Fy_modal uExtrap_modal vExtrap_modal

%% ========================================================================
%  Part 4: Plot MTU Upper Bounds vs Trajectory Errors (Figure 7)
%  ========================================================================

fprintf('\n========================================\n');
fprintf('PART 4: Plotting Results (Figure 7)\n');
fprintf('========================================\n\n');

% Load all saved data
fprintf('Loading saved data...\n');

% Load MTU bound data
d1 = load(sprintf('bounddatasave_backwards_%s_fmnm1m2_gm3.mat', pNames{1}));
d2 = load(sprintf('bounddatasave_backwards_%s_fmnm1m2_gm3.mat', pNames{2}));
d3 = load(sprintf('bounddatasave_backwards_%s_fmnm1m2_gm3.mat', pNames{3}));

% Load particle position data
p1_orig = load(sprintf('particlePos_backwards_%s_original.mat', pNames{1}));
p1_modal = load(sprintf('particlePos_backwards_%s_fmnm1m2.mat', pNames{1}));

p2_orig = load(sprintf('particlePos_backwards_%s_original.mat', pNames{2}));
p2_modal = load(sprintf('particlePos_backwards_%s_fmnm1m2.mat', pNames{2}));

p3_orig = load(sprintf('particlePos_backwards_%s_original.mat', pNames{3}));
p3_modal = load(sprintf('particlePos_backwards_%s_fmnm1m2.mat', pNames{3}));

% Extract point coordinates
px1 = p1_orig.xMinROI + p1_orig.dtol;
py1 = p1_orig.yMinROI + p1_orig.dtol;

px2 = p2_orig.xMinROI + p2_orig.dtol;
py2 = p2_orig.yMinROI + p2_orig.dtol;

px3 = p3_orig.xMinROI + p3_orig.dtol;
py3 = p3_orig.yMinROI + p3_orig.dtol;

fprintf('\nPlotting trajectory errors for points:\n');
fprintf('  Point 1: (%.3f, %.3f)\n', px1, py1);
fprintf('  Point 2: (%.3f, %.3f)\n', px2, py2);
fprintf('  Point 3: (%.3f, %.3f)\n\n', px3, py3);

% Grid indices (center point)
ii = 2;
jj = 2;

% Time parameters
tLength_plot = d1.tLength;
dt_plot = d1.dt;
tVec_plot = (0:1:abs(tLength_plot)-1) / 60;  % Normalize by oscillation period

% Compute metrics for each point
% Point 1
cseInt1 = reshape(d1.cseIntegral(ii, jj, :), [abs(tLength_plot), 1]);
delta1 = reshape(d1.deltaInfty(ii, jj, :), [abs(tLength_plot), 1]);
err1 = (p1_orig.xPos - p1_modal.xPos).^2 + (p1_orig.yPos - p1_modal.yPos).^2;

% Point 2
cseInt2 = reshape(d2.cseIntegral(ii, jj, :), [abs(tLength_plot), 1]);
delta2 = reshape(d2.deltaInfty(ii, jj, :), [abs(tLength_plot), 1]);
err2 = (p2_orig.xPos - p2_modal.xPos).^2 + (p2_orig.yPos - p2_modal.yPos).^2;

% Point 3
cseInt3 = reshape(d3.cseIntegral(ii, jj, :), [abs(tLength_plot), 1]);
delta3 = reshape(d3.deltaInfty(ii, jj, :), [abs(tLength_plot), 1]);
err3 = (p3_orig.xPos - p3_modal.xPos).^2 + (p3_orig.yPos - p3_modal.yPos).^2;

% Generate figure
fprintf('Generating Figure 7...\n');

% Define colors
c1 = [0, 0, 0, 1];                   % Black for point 1
c2 = [89, 99, 235, 255] / 255;       % Blue for point 2
c3 = [235, 99, 89, 255] / 255;       % Red for point 3

% Plot settings
SetFont = 'Times New Roman';
axislength = 18;
ep = 0.01;

% Create figure
figure('Position', [10, 10, 470, 500]);
set(gcf, 'Name', 'Figure 7: MTU Upper Bounds vs Trajectory Errors');

% Plot upper bound estimates (solid lines)
h1 = semilogy(tVec_plot, (cseInt1 .* delta1), ...
    'color', c1, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', ...
    sprintf('$(x,y) = (%.2f, %.2f)$', px1, py1));
hold on;
h2 = semilogy(tVec_plot, (cseInt2 .* delta2), ...
    'color', c2, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', ...
    sprintf('$(x,y) = (%.2f, %.2f)$', px2, py2));
h3 = semilogy(tVec_plot, (cseInt3 .* delta3), ...
    'color', c3, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', ...
    sprintf('$(x,y) = (%.2f, %.2f)$', px3, py3));

% Plot actual trajectory errors (dotted lines)
semilogy(tVec_plot, reshape(err1(ii, jj, 1:abs(tLength_plot)), [1, abs(tLength_plot)]), ...
    'color', c1, 'Linewidth', 2, 'LineStyle', ':', 'HandleVisibility', 'off');
semilogy(tVec_plot, reshape(err2(ii, jj, 1:abs(tLength_plot)), [1, abs(tLength_plot)]), ...
    'color', c2, 'Linewidth', 2, 'LineStyle', ':', 'HandleVisibility', 'off');
semilogy(tVec_plot, reshape(err3(ii, jj, 1:abs(tLength_plot)), [1, abs(tLength_plot)]), ...
    'color', c3, 'Linewidth', 2, 'LineStyle', ':', 'HandleVisibility', 'off');

% Set axis properties
xlim([0, 1]);
ylim([10^-8, 10^0]);
set(gca, 'LineWidth', 1.2, 'TickLength', [ep, ep]);
set(gcf, 'color', 'w');

% Labels
xlabel('$t/T_{p}$', 'Interpreter', 'latex', 'fontsize', axislength);
ylabel('Squared Error', 'fontsize', axislength);

% Legend
legend([h1, h2, h3], 'Location', 'southeast', 'FontSize', 18, 'Interpreter', 'latex');

% Final formatting
ax = gca;
box on;
ax.FontName = SetFont;
ax.FontSize = axislength;
hold off;

fprintf('\nFigure 7 generated successfully!\n');
fprintf('Solid lines: Upper bound estimates\n');
fprintf('Dotted lines: Actual trajectory errors\n\n');

fprintf('========================================\n');
fprintf('Complete MTU Analysis Finished!\n');
fprintf('========================================\n');
fprintf('\nSummary:\n');
fprintf('  - MTU bounds computed for %d points\n', length(pxList));
fprintf('  - Particle trajectories computed (original + modal)\n');
fprintf('  - Figure 7 generated\n\n');
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
        utemp = [uleft, uvec(:,:,i)];
        uExt(:,:,i) = [utb; utemp; utb];
    end
    [yi,xi,ti] = ndgrid(yVecExt,xVecExt, tVec);
    F = griddedInterpolant(yi, xi, ti, uExt,'makima','nearest'); 
end

function u = uExtrapolateM(x, y, t, u, F) 
    u = F(y,x,t);
end

