
%% ========================================================================
%  Plot Modal Representation and Perturbation Fields for Cylinder Wake
%  ========================================================================
%  This script generates Figure 5 from the paper:
%  (a) Modal representation: tilde{u} = umean + u1 + u2
%  (b) Perturbation: u' = u3
%  ========================================================================

clear all

% Add necessary paths
addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

% Load data
fprintf('Loading data...\n');
load('cylindervelocity_zen.mat');
load('cylinderPODmodesPair.mat');
load('cylcoords.mat');

%% ========================================================================
%  Setup Parameters
%  ========================================================================

% Frame to plot (first frame)
frame_idx = 1;

% Custom colormap
mycolormap = customcolormap(linspace(0, 1, 13), ...
    {'#a9212d', '#b8412a', '#ca6827', '#edb121', '#f5d586', ...
     '#ffffff', '#ffffff', '#ffffff', ...
     '#a9d9f8', '#2295e0', '#0758ab', '#0d449d', '#132d8d'});

% Plot parameters
width = 800;
height = 320;
axislength = 20;
ep = 0.01;
SetFont = 'Times New Roman';

% Spatial limits
xlim_range = [-2, 20];
ylim_range = [-4, 4];

%% ========================================================================
%  Figure 5(a): Modal Representation (umean + u1 + u2)
%  ========================================================================

fprintf('\nPlotting modal representation field...\n');

% Compute modal representation: tilde{u} = umean + u1 + u2
u_modal = umean + u1(:, :, frame_idx) + u2(:, :, frame_idx);

% Plot settings for modal representation
umin_modal = 0;
umax_modal = 1.2;

figure(1);
set(gcf, 'Position', [100, 100, width, height]);
contourf(xMat, yMat, u_modal, 80, 'LineStyle', 'none');
colormap(mycolormap);
axis equal;
clim([umin_modal, umax_modal]);

% Colorbar
c = colorbar;
c.Ticks = linspace(umin_modal, umax_modal, 3);
c.FontSize = axislength;
c.FontName = SetFont;
c.LineWidth = 1.2;
c.TickLength = 0.02;

hold on;

% Add cylinder
fill(xcyl, ycyl, [0.6, 0.6, 0.6]);
plot(xcyl, ycyl, 'k', 'LineWidth', 1.2);

% Labels and formatting
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', axislength);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', axislength);
title('(a) Modal Representation: $\tilde{u} = \bar{u} + u_1 + u_2$', ...
    'Interpreter', 'latex', 'fontsize', axislength);
xlim(xlim_range);
ylim(ylim_range);
xticks(0:4:20);
yticks(-4:4:4);

ax = gca;
box on;
ax.FontName = SetFont;
ax.FontSize = axislength;
ax.LineWidth = 1.2;
ax.TickLength = [ep, ep];
set(gcf, 'color', 'w');

hold off;

fprintf('Modal representation plot complete.\n');

%% ========================================================================
%  Figure 5(b): Perturbation (u3)
%  ========================================================================

fprintf('\nPlotting perturbation field...\n');

% Perturbation: u' = u3
u_pert = u3(:, :, frame_idx);

% Plot settings for perturbation
umin_pert = -0.04;
umax_pert = 0.04;

figure(2);
set(gcf, 'Position', [100, 500, width, height]);
contourf(xMat, yMat, u_pert, 100, 'LineStyle', 'none');
colormap(mycolormap);
axis equal;
clim([umin_pert, umax_pert]);

% Colorbar
c = colorbar;
c.Ticks = linspace(umin_pert, umax_pert, 3);
c.FontSize = axislength;
c.FontName = SetFont;
c.LineWidth = 1.2;
c.TickLength = 0.02;

hold on;

% Add cylinder
fill(xcyl, ycyl, [0.3, 0.3, 0.3]);
plot(xcyl, ycyl, 'k', 'LineWidth', 1.2);

% Labels and formatting
xlabel('$x$', 'Interpreter', 'latex', 'fontsize', axislength);
ylabel('$y$', 'Interpreter', 'latex', 'fontsize', axislength);
title('(b) Perturbation: $u^\prime = u_3$', ...
    'Interpreter', 'latex', 'fontsize', axislength);
xlim(xlim_range);
ylim(ylim_range);
xticks(0:4:20);
yticks(-4:4:4);

ax = gca;
box on;
ax.FontName = SetFont;
ax.FontSize = axislength;
ax.LineWidth = 1.2;
ax.TickLength = [ep, ep];
set(gcf, 'color', 'w');

hold off;

fprintf('Perturbation plot complete.\n');
fprintf('\n========================================\n');
fprintf('Both plots generated successfully!\n');
fprintf('Figure 1: Modal representation\n');
fprintf('Figure 2: Perturbation\n');
fprintf('========================================\n');