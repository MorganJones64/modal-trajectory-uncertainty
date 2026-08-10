%% DMD Recon Flow Field
% C:\Users\Gator\Documents\PhD Research 2023\Lagrangian Coherent Structures\Haller_2020\experimental data
clear all
addpath('./data')
addpath('./customcolormap/')
load('cylindervelocity_zen.mat')
load("cylinderPODmodesPair.mat")
%%
ufq = umean+u1+u2;
vfq = vmean+v1+v2;
ugq = u3;
vgq = v3;
clear ufint vfint ugint vgint u1 u2 u3 u4 umean v1 v2 v3 v3 v4 vmean xnMat ynMat
%%
ufMesh=ufq;
vfMesh=vfq;
ugMesh=ugq;
vgMesh=vgq;
%[ny,nx,nt]=size(ur);
clear ufq vfq ugq vgq

%% Compute FTLE
%yVec=y(129:385,1);
%xVec=x(1,129:end)';
tVec = (0:1:(nt-1))*dt';
yVec=y; %y(1,:)';
xVec=x;%x(1,256:769)';

Ffx = interpolantxy(ufMesh,xVec,yVec,tVec,1);
Ffy = interpolantxy(vfMesh,xVec,yVec,tVec,0);
Fgx = interpolantxy(ugMesh,xVec,yVec,tVec,0);
Fgy = interpolantxy(vgMesh,xVec,yVec,tVec,0);

ufExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Ffx);
vfExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Ffy);
ugExtrap = @(x, y, t, u) uExtrapolateM(x, y, t, u, Fgx);
vgExtrap = @(x, y, t, v) uExtrapolateM(x, y, t, v, Fgy);
clear Ffx Ffy Fgx Fgy

% Interpolate data
mycolormap = customcolormap(linspace(0,1,11),{'#a9212d','#b8412a','#ca6827','#edb121','#f5d586','#ffffff','#a9d9f8','#2295e0','#0758ab','#0d449d','#132d8d'});

%% Compute FTLE
clear sigmaS

tol  = 1;
dtol = 1e-3;

% =========================
% Define all probe points
% =========================
% pxList = [2.02, 4.17, 6.00, 6.37];
% pyList = [-0.41, -0.44, -0.589, 0.138];
% pxList = [4.14, 3.67, 4.33];
% pyList = [0.363, 0.013, 0.04];
%Farfield region
%pxList = [9.426, 13.72, 15.677];
%pyList = [-1.42, 1.14, 1.316];
%New farfield points
pxList = [14.053, 11.189];
pyList = [-2.419, 1.942];

pNames = {'p12','p13'};

ROInx = 3;
ROIny = 3;

method = 'RK4';

%% Backwards FTLE settings

tLength = -199;   % negative if backwards integration
finc    = 1;
tStep   = -1;
fstart  = 200;

clc
addpath('fxn\')

% =====================================================
% LOOP OVER ALL POINTS
% =====================================================
for p = 1:length(pxList)

    % Current point
    px = pxList(p);
    py = pyList(p);

    fprintf('\n=====================================\n');
    fprintf('Running %s\n', pNames{p});
    fprintf('px = %.5f, py = %.5f\n', px, py);
    fprintf('=====================================\n');

    % ROI around point
    xMinROI = px - dtol;
    xMaxROI = px + dtol;

    yMinROI = py - dtol;
    yMaxROI = py + dtol;

    % ==========================================
    % Preallocate arrays
    % ==========================================
    sigma_ftle  = zeros([ROInx, ROIny, abs(tLength)]);
    cseIntegral = zeros([ROInx, ROIny, abs(tLength)]);
    deltaInfty  = zeros([ROInx, ROIny, abs(tLength)]);

    % ==========================================
    % Metadata
    % ==========================================
    notes = sprintf('f=m1+m2+mean, g=m3, backward bound estimate, %s', ...
        pNames{p});

    filename = sprintf('bounddatasave_backwards_%s_fmnm1m2_gm3.mat', ...
        pNames{p});

    % ==========================================
    % Initial save
    % ==========================================
    save(filename, ...
        'sigma_ftle', ...
        'cseIntegral', ...
        'deltaInfty', ...
        'dt', ...
        'tLength', ...
        'xMaxROI', ...
        'xMinROI', ...
        'yMaxROI', ...
        'yMinROI', ...
        'fstart', ...
        'method', ...
        'tol', ...
        'notes');

    % ==========================================
    % Run computation
    % ==========================================
    tic

    [cseIntegral, deltaInfty, xPos, yPos] = MS2bound( ...
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

    toc

    % ==========================================
    % Append results
    % ==========================================
    save(filename, ...
        'cseIntegral', ...
        'deltaInfty', ...
        'xPos', ...
        'yPos', ...
        '-append');

    fprintf('Finished %s\n', pNames{p});

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

