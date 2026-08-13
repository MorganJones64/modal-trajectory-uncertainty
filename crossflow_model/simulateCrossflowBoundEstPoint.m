%% ========================================================================
%  Compute Point-Wise MTU Bounds for Crossflow Model
%  ========================================================================
%  Computes MTU upper bound estimates at three probe locations and
%  compares with actual trajectory errors. Generates Figure 1 from paper.
%  ========================================================================

clear all
close all
addpath('models');
addpath('data');

fprintf('\n========================================\n');
fprintf('Point-Wise MTU Bound Computation\n');
fprintf('========================================\n\n');

% Model Parameters
% ========================================================================
A = 0.1;            % Gyre amplitude
gamma = 0.25;       % Asymmetry parameter
U = 0.2;            % Background velocity
om = 0;             % Gyre oscillation frequency (steady)

epsn = 0.03;        % Jet amplitude
omega = 2*pi/5;     % Jet frequency (T_p = 5 s)
bt = 4;             % Jet width
at = 1.4;           % Jet center

dt = 0.025;         % Time step
T = 20;             % Integration duration
int = 'f';          % Forward integration

% Define three probe locations
tol = 1e-7;
locations = [
    0.48,  0.63;    % Location 1
    1.50,  0.75;    % Location 2  
    1.13,  0.918    % Location 3
];
locationNames = {'Location 1', 'Location 2', 'Location 3'};

fprintf('Model Parameters:\n');
fprintf('  Base flow: A=%.2f, U=%.2f\n', A, U);
fprintf('  Jet: eps=%.3f, omega=%.4f\n', epsn, omega);
fprintf('  Integration time: %.1f s\n', T);
fprintf('\nProbe Locations:\n');
for i = 1:size(locations, 1)
    fprintf('  %s: (%.3f, %.3f)\n', locationNames{i}, locations(i,1), locations(i,2));
end
fprintf('\n');

% Integration parameters
intLength = T/dt;
if strcmp(int, 'f')
    sgn = 1;
else
    sgn = -1;
end

% Storage for all locations
all_mseIntegral = cell(size(locations, 1), 1);
all_deltaInfty = cell(size(locations, 1), 1);
all_xT_pert = cell(size(locations, 1), 1);
all_yT_pert = cell(size(locations, 1), 1);
all_xT_clean = cell(size(locations, 1), 1);
all_yT_clean = cell(size(locations, 1), 1);



%% Loop Over All Locations
% ========================================================================
total_tic = tic;

for loc_idx = 1:size(locations, 1)   
    fprintf('\n========================================\n');
    fprintf('Processing %s: (%.3f, %.3f)\n', locationNames{loc_idx}, ...
        locations(loc_idx,1), locations(loc_idx,2));
    fprintf('========================================\n\n');
    
    % Define ROI around location
    xgmin = locations(loc_idx, 1) - tol;
    xgmax = locations(loc_idx, 1) + tol;
    ygmin = locations(loc_idx, 2) - tol;
    ygmax = locations(loc_idx, 2) + tol;
    
    % Create 3x3 grid
    xgrid = linspace(xgmin, xgmax, 3);
    ygrid = linspace(ygmin, ygmax, 3);
    [x0, y0] = meshgrid(xgrid, ygrid);
    
    dx = abs(xgrid(2) - xgrid(1));
    dy = dx;
    
    yIC(1,:,:) = x0';
    yIC(2,:,:) = y0';
    
    % Part 1: Compute MTU Bounds
    % ========================================================================
    fprintf('Part 1: Computing MTU bounds...\n');
    
    t0 = 0;
    tend = intLength + t0 - 1;
    mseIntegral = zeros(length(xgrid), length(ygrid), intLength);
    deltaInfty = zeros(length(xgrid), length(ygrid), intLength);
    
    tic;
    for t = t0:tend
        if mod(t, 100) == 0
            fprintf('  Progress: %.1f%%\n', 100*t/tend);
        end
        
        for s = t0:(t+1)
            if strcmp(int, 'f')
                sVec1 = t0:s;
                sVec2 = s:t;
            else
                sVec1 = flip(t0:s);
                sVec2 = flip(s:t);
            end
            
            yin = yIC;
            
            % Advect t0->s (base flow only)
            if s > t0
                for i = sVec1;
                    time = i*dt;
                    yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,0,0,0,0), ...
                        sgn*dt, sgn*time, yin);
                    yin = yout;
                end
            else
                time = t0*dt; %special case: s=t0, no advection needed
            end
            
            % Compute perturbation amplitude at s
            g0(:,:,:,s+1) = gVEC(time, yin, omega, bt, at);
            
            % Advect s->t (base flow only)
            for i = sVec2
                time = i*dt;
                yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,0,0,0,0), ...
                    sgn*dt, sgn*time, yin);
                yin = yout;
            end
            
            % Reshape
            xT = reshape(yout(1,:,:), length(xgrid), length(ygrid));
            yT = reshape(yout(2,:,:), length(xgrid), length(ygrid));
            
            % Gradients
            [dxTdx0, dxTdy0] = gradient(xT, dx, dx);
            [dyTdx0, dyTdy0] = gradient(yT, dx, dx);
            
            % Compute eigenvalues
            for i = 1:length(xgrid)
                for j = 1:length(ygrid)
                    D = [dxTdx0(i,j), dxTdy0(i,j); dyTdx0(i,j), dyTdy0(i,j)];
                    sigma_mse(i,j,s+1) = sqrt(max(eig(D'*D)));
                end
            end
        end
        
        % Compute MSE integral
        for i = 1:length(xgrid)
            for j = 1:length(ygrid)
                mseIntegral(i,j,t+1) = trapz(sigma_mse(i,j,:)*dt);
            end
        end
        
        % Compute deltaInfty
        gx0 = g0(1,:,:,:);
        gy0 = g0(2,:,:,:);
        deltaInfty(:,:,t+1) = epsn*reshape(max(sqrt(gx0.^2+gy0.^2),[],4), ...
            [length(xgrid),length(ygrid)]);
        
        clear g0 gx0 gy0 sigma_mse
    end
    
    loc_time = toc;
    fprintf('Part 1 complete! Time: %.2f s\n\n', loc_time);
    
    % Store results
    all_mseIntegral{loc_idx} = mseIntegral;
    all_deltaInfty{loc_idx} = deltaInfty;
    
    % Part 2: Compute Particle Trajectories (With Perturbation)
    fprintf('Part 2: Computing trajectories with perturbation...\n');
    
    xT_pert = zeros(length(xgrid), length(ygrid), intLength);
    yT_pert = zeros(length(xgrid), length(ygrid), intLength);
    
    yin = yIC;
    for t = 0:(intLength-1)
        time = t*dt;
        yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,epsn,omega,bt,at), ...
            sgn*dt, sgn*time, yin);
        yin = yout;
        xT_pert(:,:,t+1) = reshape(yout(1,:,:), length(xgrid), length(ygrid));
        yT_pert(:,:,t+1) = reshape(yout(2,:,:), length(xgrid), length(ygrid));
    end
    
    all_xT_pert{loc_idx} = xT_pert;
    all_yT_pert{loc_idx} = yT_pert;
    
    fprintf('Part 2 complete!\n\n');
    
    % Part 3: Compute Particle Trajectories (Without Perturbation)
    fprintf('Part 3: Computing trajectories without perturbation...\n');
    
    xT_clean = zeros(length(xgrid), length(ygrid), intLength);
    yT_clean = zeros(length(xgrid), length(ygrid), intLength);
    
    yin = yIC;
    for t = 0:(intLength-1)
        time = t*dt;
        yout = rk4singlestep(@(t,x)crossflowVEC(t,x,U,A,gamma,om,0,0,0,0), ...
            sgn*dt, sgn*time, yin);
        yin = yout;
        xT_clean(:,:,t+1) = reshape(yout(1,:,:), length(xgrid), length(ygrid));
        yT_clean(:,:,t+1) = reshape(yout(2,:,:), length(xgrid), length(ygrid));
    end
    
    all_xT_clean{loc_idx} = xT_clean;
    all_yT_clean{loc_idx} = yT_clean;
    
    fprintf('Part 3 complete!\n');
    
end

total_time = toc(total_tic);

fprintf('\n========================================\n');
fprintf('All Computations Complete!\n');
fprintf('Total time: %.2f seconds\n', total_time);
fprintf('========================================\n\n');

% Save Results
filename = 'data/crossflow_point_bounds.mat';
save(filename, 'all_mseIntegral', 'all_deltaInfty', ...
    'all_xT_pert', 'all_yT_pert', 'all_xT_clean', 'all_yT_clean', ...
    'locations', 'locationNames', 'dt', 'intLength', 'epsn', 'T');
fprintf('Results saved to: %s\n\n', filename);

%% Plot Results
fprintf('Generating plots...\n');

% Colors
c1 = [44, 153, 61, 255]/255;
c2 = [89, 99, 235, 255]/255;
c3 = [235, 99, 89, 255]/255;

SetFont = 'Times New Roman';
axislength = 18;
ep = 0.01;
tVec = (0:1:intLength-1)*dt;

ii = 2;  % Center point
jj = 2;



for loc_idx = 1:size(locations, 1)
    figure(loc_idx);
    set(gcf, 'Position', [100, 100, 420, 470]);
    % Extract data
    cseInt = reshape(all_mseIntegral{loc_idx}(ii,jj,:), [intLength,1]);
    delta = reshape(all_deltaInfty{loc_idx}(ii,jj,:), [intLength,1]);
    err = ((all_xT_pert{loc_idx} - all_xT_clean{loc_idx}).^2 + ...
           (all_yT_pert{loc_idx} - all_yT_clean{loc_idx}).^2).^2;
    
    % Select color
    if loc_idx == 1
        col = c1;
    elseif loc_idx == 2
        col = c2;
    else
        col = c3;
    end
    
    % Plot upper bound (solid line)
    semilogy(tVec, (cseInt.*delta).^2, 'Color', col, 'LineWidth', 2, ...
        'LineStyle', '-', 'DisplayName', sprintf('Upper Bound: (%.3f, %.3f)', ...
        locations(loc_idx,1), locations(loc_idx,2)));
    hold on;
    
    % Plot actual error (dotted line)
    center_error = reshape(err(ii,jj,:), [1,intLength]);
    semilogy(tVec, center_error, 'Color', col, 'LineWidth', 2, ...
        'LineStyle', ':', 'HandleVisibility', 'off');
    
    xlim([0, T]);
    ylim([10^-15, 10^5]);
    set(gca, 'LineWidth', 1.2, 'TickLength', [ep, ep]);
    set(gcf, 'color', 'w');
    xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', axislength);
    ylabel('Squared Error', 'FontSize', axislength);
    legend('Location', 'northwest', 'FontSize', 14);
    ax = gca;
    box on;
    ax.FontName = SetFont;
    ax.FontSize = axislength;
    hold off;
end

fprintf('\n========================================\n');
fprintf('Point-Wise MTU Analysis Complete!\n');
fprintf('========================================\n');
fprintf('\nFigure: Upper bounds (solid) vs actual errors (dotted)\n');
fprintf('  - Solid lines: MTU upper bound estimates\n');
fprintf('  - Dotted lines: Actual trajectory errors\n\n');