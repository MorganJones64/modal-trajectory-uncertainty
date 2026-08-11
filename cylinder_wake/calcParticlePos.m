%% ========================================================================
%  Compute Particle Positions for Cylinder Wake
%  ========================================================================
%  This script computes particle trajectories at specific probe points
%  for both original and modal representation flow fields.
%  Used to generate trajectory error plots (Figure 7 in paper).
%  ========================================================================

clear all

% Add necessary paths
addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

fprintf('\n========================================\n');
fprintf('Particle Position Computation\n');
fprintf('========================================\n\n');

%% ========================================================================
%  Load Data
%  ========================================================================

fprintf('Loading data...\n');
load('cylindervelocity_zen.mat');
load('cylinderPODmodesPair.mat');

% Use original flow field for particle tracking
uMesh = u;
vMesh = v;

%% ========================================================================
%  Define Probe Points
%  ========================================================================

% Probe point parameters
dtol = 1e-3;  % Small tolerance around each point
ROInx = 3;    % 3x3 grid
ROIny = 3;
method = 'RK4';

% Define probe point locations (should match calcMTUboundSweep.m)
pxList = [15.68, 11.189, 4.33];
pyList = [1.32, 1.942, 0.04];
pNames = {'p11', 'p13', 'p6c'};

fprintf('\nProbe points:\n');
for i = 1:length(pxList)
    fprintf('  %s: (x, y) = (%.3f, %.3f)\n', pNames{i}, pxList(i), pyList(i));
end
fprintf('\n');

%% ========================================================================
%  Setup Computation Parameters
%  ========================================================================

% Time parameters
tVec = (0:1:nt-1) * dt';
yVec = y;
xVec = x;

% Backward FTLE parameters
tLength = -199;  % Negative for backward integration
finc = 1;
tStep = -1;
fstart = 200;

fprintf('Integration parameters:\n');
fprintf('  Integration length: %d frames\n', abs(tLength));
fprintf('  Starting frame: %d\n', fstart);
fprintf('  Method: %s\n\n', method);

% Create extrapolation functions
fprintf('Creating interpolation functions...\n');
Fx = interpolantxy(uMesh, xVec, yVec, tVec, 1);
Fy = interpolantxy(vMesh, xVec, yVec, tVec, 0);

uExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx);
vExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy);

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
    
    % Output filename
    filename = sprintf('particlePos_backwards_%s_original.mat', pNames{p});
    
    % Run FTLE particle advection
    fprintf('Computing particle trajectories...\n');
    tic;
    
    [sigma_ftle, xPos, yPos] = FTLEsb( ...
        uMesh, vMesh, ...
        xVec, yVec, ...
        fstart, tLength, tStep, dt, ...
        xMinROI, xMaxROI, ...
        yMinROI, yMaxROI, ...
        ROInx, ROIny, method, ...
        'extrap', true, ...
        'uExtrap', uExtrap, ...
        'vExtrap', vExtrap);
    
    elapsed = toc;
    fprintf('Computation time: %.2f seconds\n', elapsed);
    
    % Save results
    save(filename, ...
        'xPos', 'yPos', 'sigma_ftle', ...
        'yMaxROI', 'yMinROI', 'xMaxROI', 'xMinROI', ...
        'tLength', 'tVec', 'dt', 'fstart', 'finc', 'dtol');
    
    fprintf('Saved: %s\n', filename);
    
end

fprintf('\n========================================\n');
fprintf('All particle position computations complete!\n');
fprintf('========================================\n\n');
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

