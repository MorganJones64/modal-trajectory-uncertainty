% C:\Users\Gator\Documents\PhD Research 2023\Lagrangian Coherent Structures\Haller_2020\cylinder wake
%% Get DAta
clear all
addpath('./data')
addpath('./fxn')
addpath('./customcolormap/')
load('cylindervelocity_zen.mat')
load("cylinderPODmodesPair.mat")
%%
%C:\Users\Gator\Documents\PhD Research 2023\Lagrangian Coherent Structures\Haller_2020\cylinder wakeload("cylinderDMDmodes.mat")
% Interpolate data
uq = u;
vq = v;
uMesh=uq;
vMesh=vq;
% ==========================================
% Specify regions / probe points
% ==========================================
dtol = 1e-3;

% pxList = [2.02, 4.17, 6.00, 6.37];
% pyList = [-0.41, -0.44, -0.589, 0.138];
% 
% pxList = [2.02, 4.14, 6.00, 6.37];
% pyList = [-0.41, 0.363, -0.589, 0.138];
% 
% pNames = {'p5','p6','p7','p8'};

% pxList = [4.14, 3.67, 4.33];
% pyList = [0.363, 0.013, 0.04];
% 
% pNames = {'p6a','p6b','p6c'};

% pxList = [9.426, 13.72, 15.677];
% pyList = [-1.42, 1.14, 1.316];
% 
% pNames = {'p9','p10','p11'};

pxList = [14.053, 11.189];
pyList = [-2.419, 1.942];

pNames = {'p12','p13'};

ROInx = 3;
ROIny = 3;

method = 'RK4';

% ==========================================
% Time vectors
% ==========================================
tVec = (0:1:nt-1)*dt';

yVec = y;
xVec = x;

% ==========================================
% Extrapolation functions
% ==========================================
Fx = interpolantxy(uMesh, xVec, yVec, tVec, 1);
Fy = interpolantxy(vMesh, xVec, yVec, tVec, 0);

uExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fx);
vExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fy);

clc

%% Backwards FTLE settings

tLength = -199;
finc    = 1;
tStep   = -1;
fstart  = 200;

% ==========================================
% LOOP OVER ALL POINTS
% ==========================================
for p = 1:length(pxList)

    % Current point
    px = pxList(p);
    py = pyList(p);

    fprintf('\n=====================================\n');
    fprintf('Running %s\n', pNames{p});
    fprintf('px = %.5f, py = %.5f\n', px, py);
    fprintf('=====================================\n');

    % ROI
    xMinROI = px - dtol;
    xMaxROI = px + dtol;

    yMinROI = py - dtol;
    yMaxROI = py + dtol;

    % Output filename
    filename = sprintf('particlePos_backwards_%s_original_.mat', ...
        pNames{p});

    % ==========================================
    % Run FTLE particle advection
    % ==========================================
    tic

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

    toc

    % ==========================================
    % Save results
    % ==========================================
    save(filename, ...
        'xPos', ...
        'yPos', ...
        'sigma_ftle',...
        'yMaxROI', ...
        'yMinROI', ...
        'xMaxROI', ...
        'xMinROI', ...
        'tLength', ...
        'tVec', ...
        'dt', ...
        'fstart', ...
        'finc', ...
        'dtol');

    fprintf("saved: %s\n", filename);

end
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

