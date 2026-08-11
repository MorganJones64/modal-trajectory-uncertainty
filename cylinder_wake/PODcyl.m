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

% Define truncation levels to visualize
sv1 = 2;  % First mode pair
sv2 = 4;  % First two mode pairs
sv3 = 8;  % First four mode pairs
d = 40;   % Text offset for labels

% Extract singular values
sigma = diag(Sigma);

% Create figure
figure('Position', [100, 100, 1200, 500]);

% Left subplot: Singular values (log scale)
subplot(1, 2, 1);
semilogy(sigma, 'ok', 'MarkerFaceColor', 'k');
hold on;

% Highlight specific truncation levels
semilogy(sv1, sigma(sv1), '.r', 'MarkerSize', 25);
text(sv1 + d, sigma(sv1), ['r = ', num2str(sv1)], ...
    'interpreter', 'latex', 'fontsize', 14, 'Color', 'r');

semilogy(sv2, sigma(sv2), '.r', 'MarkerSize', 25);
text(sv2 + d, sigma(sv2), ['r = ', num2str(sv2)], ...
    'interpreter', 'latex', 'fontsize', 14, 'Color', 'r');

semilogy(sv3, sigma(sv3), '.r', 'MarkerSize', 25);
text(sv3 + d, sigma(sv3), ['r = ', num2str(sv3)], ...
    'interpreter', 'latex', 'fontsize', 14, 'Color', 'r');

ylabel('\bf{Singular value}, $\sigma_k$', 'interpreter', 'latex', 'fontsize', 14);
xlabel('$k$', 'interpreter', 'latex', 'fontsize', 14);
title('Singular Values', 'interpreter', 'latex', 'fontsize', 16);
grid on;
hold off;

% Right subplot: Cumulative energy
subplot(1, 2, 2);
plot(cumsum(sigma) / sum(sigma), 'ok', 'MarkerFaceColor', 'k');
hold on;

% Compute cumulative energy for truncation levels
sc1 = cumsum(sigma(1:sv1)) / sum(sigma);
sc2 = cumsum(sigma(1:sv2)) / sum(sigma);
sc3 = cumsum(sigma(1:sv3)) / sum(sigma);

% Highlight specific truncation levels
plot(sv1, sc1(end), '.r', 'MarkerSize', 25);
text(sv1 + d, sc1(end), ['r = ', num2str(sv1)], ...
    'interpreter', 'latex', 'fontsize', 14, 'Color', 'r');

plot(sv2, sc2(end), '.r', 'MarkerSize', 25);
text(sv2 + d, sc2(end), ['r = ', num2str(sv2)], ...
    'interpreter', 'latex', 'fontsize', 14, 'Color', 'r');

plot(sv3, sc3(end), '.r', 'MarkerSize', 25);
text(sv3 + d, sc3(end), ['r = ', num2str(sv3)], ...
    'interpreter', 'latex', 'fontsize', 14, 'Color', 'r');

ylabel('\bf{Cumulative energy}', 'interpreter', 'latex', 'fontsize', 14);
xlabel('$k$', 'interpreter', 'latex', 'fontsize', 14);
title('Energy Content', 'interpreter', 'latex', 'fontsize', 16);
set(gcf, 'Color', 'w');
grid on;
hold off;

fprintf('\nEnergy captured by different truncations:\n');
fprintf('  r = %d: %.2f%%\n', sv1, sc1(end) * 100);
fprintf('  r = %d: %.2f%%\n', sv2, sc2(end) * 100);
fprintf('  r = %d: %.2f%%\n\n', sv3, sc3(end) * 100);

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