%% ========================================================================
%  Compute Backward FTLE Fields for the Oscillating Foil Wake
%  ========================================================================
%  This script computes backward FTLE fields for:
%  1. The original PIV velocity field
%  2. The opt-DMD baseline: tilde{u} = umean + u1
%  3. The perturbed modal system: tilde{u} + u' = umean + u1 + u2
%  ========================================================================

clear all

addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

fprintf('Loading oscillating foil data and opt-DMD modes...\n');
load('data/oscfoil_data.mat');
load('data/opdmd_modes.mat');
load('data/foilcoords.mat');

%% Computation parameters
xMinROI = min(x);
xMaxROI = max(x);
yMinROI = min(y);
yMaxROI = max(y);
ROInx = 200;
ROIny = 100;
method = 'RK4';
xVec = x(:)';
yVec = y(:)';

pnt = 99;
tLength = -round(0.5*pnt);
tStep = -1;
fstart = 2*pnt;

fprintf('\nFTLE parameters:\n');
fprintf('  Start frame: %d\n', fstart);
fprintf('  Integration length: %d frames (backward)\n', abs(tLength));
fprintf('  ROI grid: %d x %d\n', ROInx, ROIny);
fprintf('  Method: %s\n\n', method);

%% Part 1: original velocity field
fprintf('Part 1/3: Computing original-flow FTLE...\n');
[sigma_original, xPos_original, yPos_original] = computeFTLEField( ...
    u, v, xVec, yVec, tVec, fstart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ROInx, ROIny, method, 1, 0);

save('data/bFTLE_original.mat', 'sigma_original', 'xPos_original', ...
    'yPos_original', 'tLength', 'dt', 'fstart', 'method', ...
    'xMinROI', 'xMaxROI', 'yMinROI', 'yMaxROI');

%% Part 2: opt-DMD baseline
fprintf('Part 2/3: Computing baseline FTLE (umean + u1)...\n');
uBaseline = umean + u1;
vBaseline = vmean + v1;
[sigma_baseline, xPos_baseline, yPos_baseline] = computeFTLEField( ...
    uBaseline, vBaseline, xVec, yVec, tVec, fstart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ROInx, ROIny, method, 1, 0);

save('data/bFTLE_baseline.mat', 'sigma_baseline', 'xPos_baseline', ...
    'yPos_baseline', 'tLength', 'dt', 'fstart', 'method', ...
    'xMinROI', 'xMaxROI', 'yMinROI', 'yMaxROI');

%% Part 3: baseline plus mode-2 perturbation
fprintf('Part 3/3: Computing perturbed-system FTLE (umean + u1 + u2)...\n');
uPerturbed = uBaseline + u2;
vPerturbed = vBaseline + v2;
[sigma_perturbed, xPos_perturbed, yPos_perturbed] = computeFTLEField( ...
    uPerturbed, vPerturbed, xVec, yVec, tVec, fstart, tLength, tStep, dt, ...
    xMinROI, xMaxROI, yMinROI, yMaxROI, ROInx, ROIny, method, 1, 0);

save('data/bFTLE_perturbed.mat', 'sigma_perturbed', 'xPos_perturbed', ...
    'yPos_perturbed', 'tLength', 'dt', 'fstart', 'method', ...
    'xMinROI', 'xMaxROI', 'yMinROI', 'yMaxROI');

%% Plot fields using the publication formatting
FTLE_original = log(sigma_original) / abs(tLength*dt);
FTLE_baseline = log(sigma_baseline) / abs(tLength*dt);
FTLE_perturbed = log(sigma_perturbed) / abs(tLength*dt);

xg = xPos_original(:,:,1);
yg = yPos_original(:,:,1);
xt = 2;
yt = max(y)/2;
foilX = [xfoil(foilStartpx:foilEndpx), foilEnd, foilStart] - xt;
foilY = [yfoil(foilStartpx:foilEndpx,fstart+phase)', 0, 0] - yt;

cmap = customcolormap(linspace(0,1,6), ...
    {'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'});
plotFTLEField(11, xg-xt, yg-yt, FTLE_original, foilX, foilY, cmap, ...
    'Backward FTLE: Original Flow');
plotFTLEField(12, xg-xt, yg-yt, FTLE_baseline, foilX, foilY, cmap, ...
    'Backward FTLE: Baseline ($\bar{u}+u_1$)');
plotFTLEField(13, xg-xt, yg-yt, FTLE_perturbed, foilX, foilY, cmap, ...
    'Backward FTLE: Perturbed ($\bar{u}+u_1+u_2$)');

fprintf('\nFTLE analysis complete. Results saved in data/.\n');

%% Helper functions
function [sigma, xPos, yPos] = computeFTLEField(uMesh, vMesh, xVec, ...
    yVec, tVec, fstart, tLength, tStep, dt, xMinROI, xMaxROI, ...
    yMinROI, yMaxROI, ROInx, ROIny, method, uBC, vBC)

    Fx = interpolantxy(uMesh, xVec, yVec, tVec, uBC);
    Fy = interpolantxy(vMesh, xVec, yVec, tVec, vBC);
    uExtrap = @(x,y,t,u) evaluateInterpolant(x,y,t,u,Fx);
    vExtrap = @(x,y,t,v) evaluateInterpolant(x,y,t,v,Fy);

    [sigma, xPos, yPos] = FTLEsb(uMesh, vMesh, xVec, yVec, ...
        fstart, tLength, tStep, dt, xMinROI, xMaxROI, yMinROI, ...
        yMaxROI, ROInx, ROIny, method, 'extrap', true, ...
        'uExtrap', uExtrap, 'vExtrap', vExtrap);
end

function plotFTLEField(figureNumber, xg, yg, field, foilX, foilY, cmap, plotTitle)
    width = 600;
    height = 250;
    axislength = 20;
    ep = 0.01;
    umin = 0;
    umax = 1;
    SetFont = 'Times New Roman';

    figure(figureNumber); clf
    set(gcf,'Position',[100 100 width height],'Color','w')
    contourf(xg,yg,field,80,'LineStyle','none')
    hold on
    colormap(cmap)
    axis equal
    clim([umin,umax])
    c = colorbar;
    c.Ticks = linspace(umin,umax,3);
    c.FontSize = axislength;
    c.FontName = SetFont;
    c.LineWidth = 1.2;
    c.Color = 'k';
    c.TickLength = 0.02;
    fill(foilX,foilY,'k','EdgeColor','w')
    xlim([-1.3,max(xg(:))])
    xlabel('$x/c$','Interpreter','latex','FontSize',axislength)
    ylabel('$y/c$','Interpreter','latex','FontSize',axislength)
    title(plotTitle,'Interpreter','latex','FontSize',axislength)
    set(gca,'LineWidth',1.2,'TickLength',[ep,ep], ...
        'FontName',SetFont,'FontSize',axislength,'Color','w')
    box on
    hold off
    drawnow
end

function F = interpolantxy(uvec,xVec,yVec,tVec,bcs)
    nt = length(tVec);
    nx = length(xVec);
    ny = length(yVec);
    xVecExt = [xVec(1)-1, xVec];
    yVecExt = [yVec(1)-1, yVec, yVec(end)+1];
    if bcs == 0
        uleft = zeros(ny,1);
        utb = zeros(1,nx+1);
    else
        uleft = ones(ny,1);
        utb = ones(1,nx+1);
    end
    uExt = zeros(ny+2,nx+1,nt);
    for i = 1:nt
        utemp = [uleft, uvec(:,:,i)];
        uExt(:,:,i) = [utb; utemp; utb];
    end
    [yi,xi,ti] = ndgrid(yVecExt,xVecExt,tVec);
    F = griddedInterpolant(yi,xi,ti,uExt,'makima','nearest');
end

function value = evaluateInterpolant(x,y,t,~,F)
    value = F(y,x,t);
end
