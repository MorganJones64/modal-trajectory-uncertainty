%% ========================================================================
%  Proper Orthogonal Decomposition (POD) of Cylinder Wake Flow
%  ========================================================================
%  This script performs POD analysis on the processed velocity field data
%  and extracts the dominant spatial and temporal modes.
%  ========================================================================

clear all

% Add necessary paths
addpath('data/');
addpath('functions/');

% Load processed velocity data
fprintf('Loading velocity data...\n');
load('cylindervelocity_zen.mat');

%% ========================================================================
%  Compute POD Decomposition
%  ========================================================================

% Compute mean flow
umean = mean(u, 3);
vmean = mean(v, 3);

% Subtract mean to get fluctuating components
up = u - umean;
vp = v - vmean;

% Perform POD decomposition
fprintf('Computing POD decomposition...\n');
[Zeta, Sigma, Xi] = POD(up, vp, dx, dy, dt);
fprintf('POD decomposition complete.\n\n');

%% ========================================================================
%  Visualize Singular Values and Energy Content
%  ========================================================================

% Define mode pairs to visualize
sv1 = [1, 2];  % First mode pair
sv2 = [3, 4];  % Second mode pair
sv3 = [5, 6];  % Third mode pair
sv4 = [7, 8];  % Fourth mode pair

% Extract singular values
sigma = diag(Sigma);

% Plot settings
width = 800;
height = 450;
lw = 1.25;
ep = 0.01;
axislength = 20;
SetFont = 'Times New Roman';

% Define colors for mode pairs (in order)
c1 = '#cf6548';  % Orange-red (mode 1)
c2 = '#f0b73d';  % Yellow-orange (mode 2)
c3 = '#00b336';  % Green (mode 3)
c4 = '#1b66fa';  % Blue (mode 4)

% Create figure
figure('Position', [200, 0, width, height]);

% Left subplot: Singular values (log scale)
subplot(1, 2, 1);
semilogy(sigma, 'ok', 'LineWidth', lw);
hold on;

% Highlight mode pairs with different colors
semilogy(sv1, sigma(sv1), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c1, 'LineWidth', lw);
semilogy(sv2, sigma(sv2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c2, 'LineWidth', lw);
semilogy(sv3, sigma(sv3), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c3, 'LineWidth', lw);
semilogy(sv4, sigma(sv4), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c4, 'LineWidth', lw);

ax = gca;
box on;
ax.FontName = SetFont;
ax.FontSize = axislength;
ylabel('$\sigma_{k}$', 'interpreter', 'latex');
xlabel('$k$', 'interpreter', 'latex');
set(gca, 'LineWidth', 1.2, 'TickLength', [ep, 0], ...
    'YMinorGrid', 'off', 'XMinorGrid', 'off');
grid on;
hold off;

% Right subplot: Cumulative energy
subplot(1, 2, 2);
csum = cumsum(sigma) / sum(sigma);
plot(csum, 'ok', 'LineWidth', lw);
hold on;

% Highlight mode pairs
plot(sv1, csum(sv1), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c1, 'LineWidth', lw);
plot(sv2, csum(sv2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c2, 'LineWidth', lw);
plot(sv3, csum(sv3), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c3, 'LineWidth', lw);
plot(sv4, csum(sv4), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', c4, 'LineWidth', lw);

ax = gca;
box on;
ax.FontName = SetFont;
ax.FontSize = axislength;
xlabel('$k$', 'interpreter', 'latex');
ylabel('cumulative energy', 'interpreter', 'latex');
ylim([0.3, 1]);
yticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]);
set(gcf, 'Color', 'w');
set(gca, 'LineWidth', 1.2, 'TickLength', [ep, 0]);
grid on;
hold off;

fprintf('\nEnergy captured by different mode pairs:\n');
fprintf('  Modes 1-2: %.2f%%\n', csum(sv1(end)) * 100);
fprintf('  Modes 1-4: %.2f%%\n', csum(sv2(end)) * 100);
fprintf('  Modes 1-6: %.2f%%\n', csum(sv3(end)) * 100);
fprintf('  Modes 1-8: %.2f%%\n\n', csum(sv4(end)) * 100);

%% ========================================================================
%  Reconstruct Velocity Fields from POD Modes
%  ========================================================================

% Define index ranges for u and v components
nu = 1:(nx*ny);                % u-velocity indices
nv = (nx*ny)+1:2*(nx*ny);      % v-velocity indices
num = 9:nt;                     % Truncation at mode 9 and beyond

% Reconstruct velocity fields using different mode truncations
fprintf('Reconstructing velocity fields from POD modes...\n');

% Mode pairs (each pair captures a traveling wave pattern)
U1 = Zeta(:, 1:2) * Sigma(1:2, 1:2) * Xi(:, 1:2)';     % First mode pair
U2 = Zeta(:, 3:4) * Sigma(3:4, 3:4) * Xi(:, 3:4)';     % Second mode pair
U3 = Zeta(:, 5:6) * Sigma(5:6, 5:6) * Xi(:, 5:6)';     % Third mode pair
U4 = Zeta(:, 7:8) * Sigma(7:8, 7:8) * Xi(:, 7:8)';     % Fourth mode pair
Ut = Zeta(:, num) * Sigma(num, num) * Xi(:, num)';     % Higher modes (truncated)

% Reshape mode pairs back to spatial grids
u1 = reshape(U1(nu, :), [ny, nx, nt]);
u2 = reshape(U2(nu, :), [ny, nx, nt]);
u3 = reshape(U3(nu, :), [ny, nx, nt]);
u4 = reshape(U4(nu, :), [ny, nx, nt]);

v1 = reshape(U1(nv, :), [ny, nx, nt]);
v2 = reshape(U2(nv, :), [ny, nx, nt]);
v3 = reshape(U3(nv, :), [ny, nx, nt]);
v4 = reshape(U4(nv, :), [ny, nx, nt]);

ut = reshape(Ut(nu, :), [ny, nx, nt]);
vt = reshape(Ut(nv, :), [ny, nx, nt]);

% Save mode pair reconstructions
fprintf('Saving mode pair reconstructions...\n');
save('data/cylinderPODmodesPair.mat', ...
    'u1', 'u2', 'u3', 'u4', 'ut', ...
    'v1', 'v2', 'v3', 'v4', 'vt', ...
    'umean', 'vmean');

fprintf('\n========================================\n');
fprintf('POD analysis complete!\n');
fprintf('========================================\n');
fprintf('\nOutput files:\n');
fprintf('  - cylinderPODmodesPair.mat (mode pairs)\n\n');