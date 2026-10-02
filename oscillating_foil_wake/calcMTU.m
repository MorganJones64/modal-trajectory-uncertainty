%% ========================================================================
%  Compute Modal-Trajectory Uncertainty for the Oscillating Foil Wake
%  ========================================================================
%  Paper system (equation osc-sys1):
%    Baseline:     tilde{u} = umean + u1
%    Perturbation: u' = u2
%  The calculation is performed backward in time at one initial phase.
%  ========================================================================

clear all

addpath('./data/');
addpath('./functions/');
addpath('./customcolormap/');

fprintf('Loading oscillating foil data and opt-DMD modes...\n');
load('data/oscfoil_data.mat');
load('data/opdmd_modes.mat');
load('data/foilcoords.mat');

%% Construct modal fields
fprintf('Constructing modal system:\n');
fprintf('  Baseline: umean + u1\n');
fprintf('  Perturbation: u2\n');
ufMesh = umean + u1;
vfMesh = vmean + v1;
ugMesh = u2;
vgMesh = v2;

%% MTU parameters
xMinROI = min(x);
xMaxROI = max(x);
yMinROI = min(y);
yMaxROI = max(y);
ROInx = 300;
ROIny = 200;
method = 'RK4';
xVec = x(:)';
yVec = y(:)';

pnt = 99;
tLength = -round(1.5*pnt);
tStep = -1;
fstart = 2*pnt;

fprintf('\nMTU parameters:\n');
fprintf('  Start frame: %d\n', fstart);
fprintf('  Integration length: %d frames (backward)\n', abs(tLength));
fprintf('  ROI grid: %d x %d\n', ROInx, ROIny);
fprintf('  Method: %s\n\n', method);

%% Extrapolation functions
fprintf('Creating velocity interpolants...\n');
Ffx = interpolantxy(ufMesh,xVec,yVec,tVec,1);
Ffy = interpolantxy(vfMesh,xVec,yVec,tVec,0);
Fgx = interpolantxy(ugMesh,xVec,yVec,tVec,0);
Fgy = interpolantxy(vgMesh,xVec,yVec,tVec,0);

ufExtrap = @(x,y,t,u) evaluateInterpolant(x,y,t,u,Ffx);
vfExtrap = @(x,y,t,v) evaluateInterpolant(x,y,t,v,Ffy);
ugExtrap = @(x,y,t,u) evaluateInterpolant(x,y,t,u,Fgx);
vgExtrap = @(x,y,t,v) evaluateInterpolant(x,y,t,v,Fgy);

%% Compute MTU field (long computation)
fprintf('Computing backward MTU field. This may take considerable time...\n');
tic
[sigma_ftle,cseIntegral,deltaInfty,xPos,yPos] = CSEsb( ...
    ufMesh,vfMesh,ugMesh,vgMesh,xVec,yVec,fstart,tLength,tStep, ...
    fstart,eps,dt,xMinROI,xMaxROI,yMinROI,yMaxROI,ROInx,ROIny, ...
    method,'extrap',true,'ufExtrap',ufExtrap,'vfExtrap',vfExtrap, ...
    'ugExtrap',ugExtrap,'vgExtrap',vgExtrap);
elapsedTime = toc;

notes = 'Baseline: umean+u1; perturbation: u2; backward MTU';
filename = 'data/MTU_oscillating_foil.mat';
fprintf('Saving results to %s...\n',filename);
save(filename,'sigma_ftle','cseIntegral','deltaInfty','xPos','yPos', ...
    'dt','tLength','fstart','method','notes','xMinROI','xMaxROI', ...
    'yMinROI','yMaxROI','elapsedTime','-v7.3');

%% Compute and plot derived fields
FTLE = log(sigma_ftle)/abs(tLength*dt);
MTU = (deltaInfty.*cseIntegral).^2;
scaledMTU = log(MTU)/abs(2*tLength*dt);
zeta = (log(deltaInfty)+log(cseIntegral)-log(sigma_ftle)) ...
    / abs(tLength*dt);

xg = xPos(:,:,1);
yg = yPos(:,:,1);
xt = 2;
yt = max(y)/2;
foilX = [xfoil(foilStartpx:foilEndpx),foilEnd,foilStart] - xt;
foilY = [yfoil(foilStartpx:foilEndpx,fstart+phase)',0,0] - yt;
cmap = customcolormap(linspace(0,1,6), ...
    {'#54ff82','#1fe9ff','#1FB8FF','#184ADC','#000000','#000000'});

plotField(12,xg-xt,yg-yt,FTLE,foilX,foilY,cmap,[0,1], ...
    'Backward FTLE: Baseline ($\bar{u}+u_1$)');
plotField(13,xg-xt,yg-yt,zeta,foilX,foilY,cmap,[-1.2,-0.8], ...
    'FTLE Perturbation');
plotField(14,xg-xt,yg-yt,scaledMTU,foilX,foilY,cmap,[-1,0], ...
    'Scaled Modal-Trajectory Uncertainty');

fprintf('\nMTU analysis complete in %.1f seconds.\n',elapsedTime);

%% Helper functions
function plotField(figureNumber,xg,yg,field,foilX,foilY,cmap,limits,plotTitle)
    width = 600;
    height = 250;
    axislength = 20;
    ep = 0.01;
    SetFont = 'Times New Roman';

    figure(figureNumber); clf
    set(gcf,'Position',[100 100 width height],'Color','w')
    contourf(xg,yg,field,80,'LineStyle','none')
    hold on
    colormap(cmap)
    axis equal
    clim(limits)
    c = colorbar;
    c.Ticks = linspace(limits(1),limits(2),3);
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
    xVecExt = [xVec(1)-1,xVec];
    yVecExt = [yVec(1)-1,yVec,yVec(end)+1];
    if bcs == 0
        uleft = zeros(ny,1);
        utb = zeros(1,nx+1);
    else
        uleft = ones(ny,1);
        utb = ones(1,nx+1);
    end
    uExt = zeros(ny+2,nx+1,nt);
    for i = 1:nt
        utemp = [uleft,uvec(:,:,i)];
        uExt(:,:,i) = [utb;utemp;utb];
    end
    [yi,xi,ti] = ndgrid(yVecExt,xVecExt,tVec);
    F = griddedInterpolant(yi,xi,ti,uExt,'makima','nearest');
end

function value = evaluateInterpolant(x,y,t,~,F)
    value = F(y,x,t);
end
