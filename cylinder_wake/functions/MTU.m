function [sigma_ftle,cseIntegral, deltaInfty, xPos, yPos] = MTU(ufMesh, vfMesh, ugMesh, vgMesh, xVec, yVec, ...
    t0, tLength, tStep, fStart, eps, dt, xMinROI, xMaxROI, ...
    yMinROI, yMaxROI, nx, ny, method, options)
% Function for computing 2-D finite time Lypanov exponents from a
% time series of vector fields. At the start time the flow map is
% intialized and its deformation is computed by integrating the vector
% field time series. The FTLEs are determined based on the
% stretching of the of the flow map.
%(The flow map can be thought of as seeding the flow with particles at
% the start time and observing their trajectory as time passes)

% Required Inputs
% xVec: Vector of x grid values for vector field [n x 1]
% yVec: Vector of y grid values for vector field [m x 1]
% ufMesh: Matrix of x-component of velocity [m x n x p] Meshgrid style
% vfMesh: Matrix of y-component of velocity [m x n x p] Meshgrid style
% t_start: starting time step [scalar] (index)
% t_length: length of integration in time steps [scalar] (index length)
% NOTE: Must be negative for backward integration and less than or equal to
% p dimension
% t_step: increment of time steps (index steps) NOTE: Must be
% negative for backward integration
% dt: time between time steps (Seconds)
% xMinROI: Minimum x value for flow map
% xMaxROI: Maximum x value for flow map
% yMinROI: Minimum y value for flow map
% yMaxROI: Maximum y value for flow map
% nx: Number of grid points in ROI x-direction
% ny: Number of grid points in ROI y-direction
% Method: Integration method for determining particle trajectories
% (see trajectory function)

% Optional inputs
% xMask: Vector of x locations for closed polygon mask
% yMask: Vector of y locations for closed polygon mask
% extrap: Logical indicating whether values outside vector domain
%   should be extrapolated
% ufExtrap: Function handle for extrapolation function for u velocities.
%    Function should be in the form u = u(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest u velocity value
% vfExtrap: Function handle for extrapolation function for v velocities.
%    Function should be in the form v = v(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest v velocity value

% Outputs
% sigma: Field of FTLE values. [nx x ny x t_length]
% xPos: X position of particles as they are advected.
% To plot FTLE values use first time index [nx x ny x t_length]
% yPos: Y position of particles as they are advected.
% To plot FTLE values use first time index [nx x ny x t_length]

% Authors: Chase Klewicki (main), Oliver Kahn (RKF45)


arguments
    ufMesh double
    vfMesh double
    ugMesh double
    vgMesh double
    xVec double{mustBeVector}
    yVec double{mustBeVector}
    t0(1, 1) double
    tLength(1, 1) double
    tStep(1, 1) double
    fStart(1,1) double
    eps(1, 1) double
    dt(1, 1) double
    xMinROI(1, 1) double = min(xVec)
    xMaxROI(1, 1) double = max(xVec)
    yMinROI(1, 1) double = max(yVec)
    yMaxROI(1, 1) double = max(yVec)
    nx(1, 1) double{mustBeInteger} = 100
    ny(1, 1) double{mustBeInteger} = 100
    method(1, :) char ...
        {mustBeMember(method, {'Euler', 'RKF45', 'ODE45','RK4'})} = 'Euler'
    options.xMask double{mustBeVector}
    options.yMask double{mustBeVector}
    options.extrap(1, 1) {mustBeNumericOrLogical} = false
    options.ufExtrap function_handle
    options.vfExtrap function_handle
    options.ugExtrap function_handle
    options.vgExtrap function_handle
end

% If extrapolation is true
if options.extrap

    % Define the initial region of the fluid to track
    sxcoor = linspace(xMinROI, xMaxROI, nx);
    sycoor = linspace(yMinROI, yMaxROI, ny);

    % Initial position
    [xPos, yPos] = ndgrid(sxcoor, sycoor);

    % Check for mask
    if isfield(options, 'xMask') && isfield(options, 'yMask')
        % find positions inside mask and set to nan
        [in, ~] = inpolygon(xPos, yPos, options.xMask, options.yMask);
        xPos(in) = nan;
        yPos(in) = nan;
    end

    % Integration time length
    tSpan = abs(tLength) * dt;
    

    t = tLength + t0; %end interval of integration (t)
    % Create time vector for loop (t0)
    tLoop = t0:tStep:t;

    % Initalize position matrix
    xPos = repmat(xPos, 1, 1, length(tLoop));
    yPos = repmat(yPos, 1, 1, length(tLoop));

    % Define time vector (Problem) NOTE: Generic time vector. Does not need
    % to match velocity time instances. Usage is for if user has 2 or more
    % time-steps + forwards/backwards integration
    tVec = (0:sign(tStep):tStep) * dt;
    sigma_ftle = zeros(nx, ny);
    sgn = sign(tStep);
    dt = sign(tStep)*dt;
    for s = t0:tStep:t %sweeping parameter (s)
        sVec1 = t0:tStep:s;
        sVec2 = (s+1):tStep:t;
        s
        %pd=round(abs(abs(s-t0)/abs((t-t0))) * 100,1);
         %if mod(pd,1)==0
         %    pd
         %end
        %Advect from t0->s
        for i=sVec1
            % Determine the index
            tIndex = find(sVec1 == i);
            % Compute trajectories FIX tVec
            [xPos(:, :, tIndex+1), yPos(:, :, tIndex+1)] = ...
                trajectory(xVec, yVec, i, dt, ...
                ufMesh(:, :, i), ...
                vfMesh(:, :, i), ...
                xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
                options.extrap, @options.ufExtrap, @options.vfExtrap);
        end
        
        %ufMesh(:, :, i)
        % save g value at time s sVec1(end)
        [gx0(:,:,abs(s-t0)+1), gy0(:,:,abs(s-t0)+1)] = ...
         perturbation(xVec,yVec,i,...
         ugMesh(:, :, i),...
         vgMesh(:, :, i),...
         xPos(:, :, tIndex), yPos(:, :, tIndex),...
         options.extrap, @options.ugExtrap, @options.vgExtrap);

        %Advect from s->t
        for i=sVec2
            % Determine the index
            tIndex = find(sVec2 == i);
            % Compute trajectories
            [xPos(:, :, tIndex+1), yPos(:, :, tIndex+1)] = ...
                trajectory(xVec, yVec, i, dt, ...
                ufMesh(:, :, i), ...
                vfMesh(:, :, i), ...
                xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
                options.extrap, @options.ufExtrap, @options.vfExtrap);
        end

        % Determine components of Jacobian
        [dPhiXdY, dPhiXdX] = gradient(xPos(:, :, end), sycoor, sxcoor);
    
        [dPhiYdY, dPhiYdX] = gradient(yPos(:, :, end), sycoor, sxcoor);
    
        % Construct Jacobian for each element
        A1 = cat(3, dPhiXdX, dPhiXdY);
    
        A2 = cat(3, dPhiYdX, dPhiYdY);
    
        % Create 4-D matrix of jacobians
        B = cat(4, A1, A2);
    
        % Change order of dimensions for page transpose
        B = permute(B, [3, 4, 1, 2]);
    
        % Compute stretching. note that delta is a 2x2 matrix in each of the nx
        % by ny points
        delta = pagemtimes(B, pagetranspose(B));
        
        % Initialize FTLE value
        sigma_cse = zeros(abs(tLength),nx, ny);
    
        % Compute FTLE for each point
        for j = 1:numel(sigma_cse(1,:))
            sigma_cse(abs(s-t0)+1,j) = sqrt(max(eig(delta(:, :, j))));
            if s==t0
                sigma_ftle(j) = sqrt(max(eig(delta(:, :, j))));
            end
        end
        %toc
    end

    % If extrapolation is false stop integrating trajectories when
    % particle exits vector field domain
else
    % Find bounds of vector field data
    xMin = min(xVec);
    xMax = max(xVec);
    yMin = min(yVec);
    yMax = max(yVec);
    
    % Define the initial region of the fluid to track
    sxcoor = linspace(xMinROI, xMaxROI, nx);
    sycoor = linspace(yMinROI, yMaxROI, ny);
    
    % Initial position
    [xPosTemp, yPosTemp] = ndgrid(sxcoor, sycoor);
    
    % Check for mask
    if isfield(options, 'xMask') && isfield(options, 'yMask')
        % find positions inside mask and set to nan
        [in, ~] = inpolygon(xPosTemp, yPosTemp, ...
            options.xMask, options.yMask);
        xPosTemp(in) = nan;
        yPosTemp(in) = nan;
    else
        % If no mask exists set all points to be outside mask
        in = zeros(size(xPosTemp));
    end
    
    % Integration time length
    tSpan = abs(tLength) * dt;
    
    % Integration endpoints
    t = tLength + t0;
    tLoop = t0:tStep:t;
    
    % Initialize flags and storage
    calcMTU = ones(nx, ny, 'logical');
    outDomain = zeros(nx, ny, 'logical');
    
    % Initialize MTU arrays
    sigma_ftle = zeros(nx, ny);
    sigma_cse_store = zeros(nx, ny, abs(tLength));
    gx0_all = zeros(nx, ny, abs(tLength));
    gy0_all = zeros(nx, ny, abs(tLength));
    
    % Initalize position matrix
    xPos = repmat(zeros(nx, ny), 1, 1, length(tLoop));
    yPos = repmat(zeros(nx, ny), 1, 1, length(tLoop));
    
    % Set initial positions
    xPos(:, :, 1) = xPosTemp;
    yPos(:, :, 1) = yPosTemp;
    
    % Define time vector
    tVec = (0:sign(tStep):tStep) * dt;
    sgn = sign(tStep);
    dt_signed = sign(tStep) * dt;
    
    % Loop through sweeping parameter s
    for s = t0:tStep:t
        s
        
        sVec1 = t0:tStep:s;
        sVec2 = (s+1):tStep:t;
        
        % Advect from t0 to s
        xPosTemp_s = xPosTemp;
        yPosTemp_s = yPosTemp;
        
        for i = sVec1
            % Determine the index
            tIndex = find(sVec1 == i);
            
            % Compute trajectories only for particles in domain
            [xPosTemp_s(~outDomain), yPosTemp_s(~outDomain)] = ...
                trajectory(xVec, yVec, i, dt_signed, ...
                ufMesh(:, :, i), ...
                vfMesh(:, :, i), ...
                xPosTemp_s(~outDomain), yPosTemp_s(~outDomain), method, ...
                false, @NOP, @NOP);
            
            % Check for mask
            if isfield(options, 'xMask') && isfield(options, 'yMask')
                [in, ~] = inpolygon(xPosTemp_s, yPosTemp_s, ...
                    options.xMask, options.yMask);
            end
            
            % Determine particles that have exited domain
            index = (xPosTemp_s - xMin) .* (xPosTemp_s - xMax) >= 0 | ...
                (yPosTemp_s - yMin) .* (yPosTemp_s - yMax) >= 0 | in;
            
            % Update outDomain flag
            outDomain(index) = true;
        end
        
        % Save perturbation value at time s for particles still in domain
        [gx_temp, gy_temp] = perturbation(xVec, yVec, i, ...
            ugMesh(:, :, i), ...
            vgMesh(:, :, i), ...
            xPosTemp_s, yPosTemp_s, ...
            false, @NOP, @NOP);
        
        % Store perturbation values
        gx0_all(:, :, abs(s-t0)+1) = gx_temp;
        gy0_all(:, :, abs(s-t0)+1) = gy_temp;
        
        % Continue advecting from s to t
        for i = sVec2
            % Determine the index
            tIndex = find(sVec2 == i);
            
            % Compute trajectories only for particles in domain
            [xPosTemp_s(~outDomain), yPosTemp_s(~outDomain)] = ...
                trajectory(xVec, yVec, i, dt_signed, ...
                ufMesh(:, :, i), ...
                vfMesh(:, :, i), ...
                xPosTemp_s(~outDomain), yPosTemp_s(~outDomain), method, ...
                false, @NOP, @NOP);
            
            % Check for mask
            if isfield(options, 'xMask') && isfield(options, 'yMask')
                [in, ~] = inpolygon(xPosTemp_s, yPosTemp_s, ...
                    options.xMask, options.yMask);
            end
            
            % Determine particles that have exited domain
            index = (xPosTemp_s - xMin) .* (xPosTemp_s - xMax) >= 0 | ...
                (yPosTemp_s - yMin) .* (yPosTemp_s - yMax) >= 0 | in;
            
            % Update outDomain flag
            outDomain(index) = true;
        end
        
        % Find adjacent points that need MTU calculation
        M = zeros(size(xPosTemp_s));
        M(outDomain) = 1;
        index = conv2(M, [1, 1, 1; 1, 1, 1; 1, 1, 1], 'same') > 0;
        
        % Find row and column of points
        [ix, iy] = find(index);
        
        % Compute stretching for points that exited or are adjacent
        for i = 1:numel(ix)
            if calcMTU(ix(i), iy(i))
                % Compute Jacobian components
                if (ix(i) - 1) * (ix(i) - nx) < 0 && (iy(i) - 1) * (iy(i) - ny) < 0
                    dPhiXdX = (xPosTemp_s(ix(i)+1, iy(i)) - xPosTemp_s(ix(i)-1, iy(i))) ...
                        / (sxcoor(ix(i)+1) - sxcoor(ix(i)-1));
                    dPhiXdY = (xPosTemp_s(ix(i), iy(i)+1) - xPosTemp_s(ix(i), iy(i)-1)) ...
                        / (sycoor(iy(i)+1) - sycoor(iy(i)-1));
                    dPhiYdX = (yPosTemp_s(ix(i)+1, iy(i)) - yPosTemp_s(ix(i)-1, iy(i))) ...
                        / (sxcoor(ix(i)+1) - sxcoor(ix(i)-1));
                    dPhiYdY = (yPosTemp_s(ix(i), iy(i)+1) - yPosTemp_s(ix(i), iy(i)-1)) ...
                        / (sycoor(iy(i)+1) - sycoor(iy(i)-1));
                    
                    % Construct Jacobian
                    A = [dPhiXdX, dPhiXdY; dPhiYdX, dPhiYdY];
                    delta = A' * A;
                    
                    % Compute stretching
                    sigma_cse_store(ix(i), iy(i), abs(s-t0)+1) = sqrt(max(eig(delta)));
                    
                    % Store FTLE at initial time
                    if s == t0
                        sigma_ftle(ix(i), iy(i)) = sqrt(max(eig(delta)));
                    end
                end
            end
        end
        
        % Mark these points as calculated
        calcMTU(index) = false;
    end
    
    % Calculate MTU for remaining points that stayed in domain
    [ix, iy] = find(calcMTU);
    
    % Final advection for remaining particles
    xPosTemp_final = xPosTemp;
    yPosTemp_final = yPosTemp;
    
    for i = t0:tStep:t
        [xPosTemp_final(calcMTU), yPosTemp_final(calcMTU)] = ...
            trajectory(xVec, yVec, i, dt_signed, ...
            ufMesh(:, :, i), ...
            vfMesh(:, :, i), ...
            xPosTemp_final(calcMTU), yPosTemp_final(calcMTU), method, ...
            false, @NOP, @NOP);
    end
    
    % Compute stretching for remaining points
    for i = 1:numel(ix)
        if (ix(i) - 1) * (ix(i) - nx) < 0 && (iy(i) - 1) * (iy(i) - ny) < 0
            dPhiXdX = (xPosTemp_final(ix(i)+1, iy(i)) - xPosTemp_final(ix(i)-1, iy(i))) ...
                / (sxcoor(ix(i)+1) - sxcoor(ix(i)-1));
            dPhiXdY = (xPosTemp_final(ix(i), iy(i)+1) - xPosTemp_final(ix(i), iy(i)-1)) ...
                / (sycoor(iy(i)+1) - sycoor(iy(i)-1));
            dPhiYdX = (yPosTemp_final(ix(i)+1, iy(i)) - yPosTemp_final(ix(i)-1, iy(i))) ...
                / (sxcoor(ix(i)+1) - sxcoor(ix(i)-1));
            dPhiYdY = (yPosTemp_final(ix(i), iy(i)+1) - yPosTemp_final(ix(i), iy(i)-1)) ...
                / (sycoor(iy(i)+1) - sycoor(iy(i)-1));
            
            % Construct Jacobian
            A = [dPhiXdX, dPhiXdY; dPhiYdX, dPhiYdY];
            delta = A' * A;
            
            % Store final stretching
            sigma_cse_store(ix(i), iy(i), :) = sqrt(max(eig(delta)));
            sigma_ftle(ix(i), iy(i)) = sqrt(max(eig(delta)));
        end
    end
    
    % Use stored values for output
    sigma_cse = sigma_cse_store;
    gx0 = gx0_all;
    gy0 = gy0_all;
    xPos = xPosTemp_final;
    yPos = yPosTemp_final;
    
end

cseIntegral = zeros(nx, ny);
deltaInfty = zeros(nx, ny);

%for the specific end interval of integration t, compute the CSE across the whole grid
for j = 1:numel(sigma_cse(1,:))
    cseIntegral(j) = trapz(sigma_cse(:,j)*abs(dt));
end

%g0(1,:,:,:)=gx0;
%g0(2,:,:,:)=gy0;

deltaInfty(:,:) = max(sqrt(gx0.^2 + gy0.^2),[],3); %compute amplitude of mse
%deltaInfty(:,:) = reshape(max(max(abs(g0),[],4),[],1),[nx,ny]); %compute amplitude of mse
ss=t0-fStart+1
end

function [X, Y] = trajectory(xVec, yVec, tVec, dt, uMesh, vMesh, x0, y0, ...
    method, extrapolate, ufExtrap, vfExtrap)
% Function which computes trajectory of particles in a grid by
% integrating vector field of velocity values.

% Inputs
% xVec: Vector of x grid values for vector field [n x 1] (Must be in
%   ascending order)
% yVec: Vector of y grid values for vector field [m x 1] (Must be in
%   ascending order)
% tVec: Vector of time values for vector field  [p x 1]
%   (End time can be negative or positive.)
% uMesh: Matrix of x-component of velocity [m x n x p] Meshgrid style
% vMesh: Matrix of y-component of velocity [m x n x p] Meshgrid style
% x0: Matrix of intial x-position [m x n]
% y0: Matrix of intial y-position [m x n]
% method: 
%   'Euler': A forward Euler integration scheme using cubic
%   interpolation of 2-D velocity field and tVec as integration interval
%   (~5x Faster than 'RKF45' with error O(h^2))
%   'RKF45': Runge–Kutta–Fehlberg method, error O(h^5)
%   'ODE45': Uses ode45 Matlab function and 3-D cubic interpolation to
%   compute 4th order Runge Kutta integration.
%   (extremely slow, error O(h^5))
% extrapolate: Logical indicating whether values outside vector domain
%   should be extrapolated
% ufExtrap: Function handle for extrapolation function for u velocities.
%    Function should be in the form u = u(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest u velocity value
% vfExtrap: Function handle for extrapolation function for v velocities.
%    Function should be in the form v = v(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest v velocity value

%
% Outputs
% X: x-component of trajectory [m x n]
% Y: y-component of trajectory [m x n]

arguments
    xVec double{mustBeVector}
    yVec double{mustBeVector}
    tVec double
    dt double
    uMesh double
    vMesh double
    x0 double
    y0 double
    method(1, :) char ...
        {mustBeMember(method, {'Euler', 'RKF45', 'ODE45','RK4'})} = 'Euler'
    extrapolate(1, 1) {mustBeNumericOrLogical} = false
    ufExtrap function_handle = @NOP
    vfExtrap function_handle = @NOP
end

% Initialize length
X = zeros(size(x0));
Y = zeros(size(y0));


% NaNs propogate NanNs due to cubic interpolation! Highly suggested
% that all NaNs are removed from data prior to running this function.
% Replace NaNs with mean velocity.
uMesh(isnan(uMesh)) = 0;
vMesh(isnan(vMesh)) = 0;

%sort order for forwards/backwards integration 
%[tVecSorted, indices] = sort(tVec);

%Define interpolation functions
uInterp = griddedInterpolant({xVec, yVec}, ...
    permute(uMesh(:, :), [2, 1]), ...
    'linear', 'nearest'); %changed to linear to prevent spam
vInterp = griddedInterpolant({xVec, yVec}, ...
    permute(vMesh(:, :), [2, 1]), ...
    'linear', 'nearest');


    function u = uFunction(t, x, y)

        % Interpolate to velocity find values
        % (nearest value if outside domain)
        %u = uInterp(x, y, t);
        u = uInterp(x, y);

        % Create polygon of domain
        domain = [xVec(1), yVec(1); ...
            xVec(1), yVec(end); ...
            xVec(end), yVec(end); ...
            xVec(end), yVec(1); ...
            xVec(1), yVec(1)];

        % Determine values inside domain
        [in, on] = inpolygon(x, y, domain(:, 1), domain(:, 2));

        % Conditions to be outside domain
        conditions = ~in & ~on & ~isnan(x);

        % If extrapolation is true, make the velocity values the mean
        if extrapolate
            u(conditions) = ufExtrap(x(conditions), y(conditions), ...
                t(conditions), u(conditions));
            % Else do not advect particles outside vector field domain
        else
            u(conditions) = 0;
        end

    end

    function v = vFunction(t, x, y)
        % Interpolate to velocity find values
        % (nearest value if outside domain)
        %v = vInterp(x, y, t);
         v = vInterp(x, y);

        % Create polygon of domain
        domain = [xVec(1), yVec(1); ...
            xVec(1), yVec(end); ...
            xVec(end), yVec(end); ...
            xVec(end), yVec(1); ...
            xVec(1), yVec(1)];

        % Determine values inside domain
        [in, on] = inpolygon(x, y, domain(:, 1), domain(:, 2));

        % Conditions to be in domain
        conditions = ~in & ~on & ~isnan(x);

        % If extrapolation is true
        if extrapolate
            v(conditions) = vfExtrap(x(conditions), y(conditions), ...
                t(conditions), v(conditions));
            % Else do not advect particle
        else
            v(conditions) = 0;
        end

    end


% Compute trajectory using Euler method with integration step
% equal to time vector and 2-D cubic interpolation
if strcmp(method, 'Euler')
    for i = 1:length(tVec)
        x0 = x0 + ...
            tVec(i) * uFunction(tVec(i)*ones(size(x0)), x0, y0);
        y0 = y0 + ...
            tVec(i) * vFunction(tVec(i)*ones(size(x0)), x0, y0);
    end
    X = x0;
    Y = y0;
end

%RK4 Implementation
if strcmp(method, 'RK4')
    sgn = sign(dt);
    dt = abs(dt); %abs(tVec(2)-tVec(1));
    t=tVec;

    k1x= uFunction(((t*dt)*ones(size(x0))), x0, y0);
    k1y= vFunction(((t*dt)*ones(size(x0))), x0, y0);

    k2x= uFunction(((t*dt)*ones(size(x0))+dt/2), x0 + dt/2*k1x, y0 + dt/2*k1x);
    k2y= vFunction(((t*dt)*ones(size(x0))+dt/2), x0 + dt/2*k1y, y0 + dt/2*k1y);
    
    k3x= uFunction(((t*dt)*ones(size(x0))+dt/2), x0 + dt/2*k2x, y0 + dt/2*k2x);
    k3y= vFunction(((t*dt)*ones(size(x0))+dt/2), x0 + dt/2*k2y, y0 + dt/2*k2y);

    k4x= uFunction(((t*dt)*ones(size(x0))+dt), x0 + dt*k3x, y0 + dt*k3x);
    k4y= vFunction(((t*dt)*ones(size(x0))+dt), x0 + dt*k3y, y0 + dt*k3y);

    x0 = x0 + sgn*(dt/6)*(k1x+2*k2x+2*k3x+k4x);
    y0 = y0 + sgn*(dt/6)*(k1y+2*k2y+2*k3y+k4y);
    X = x0;
    Y = y0;
end

% Compute trajectory using ODE45 and 3-D cubic interpolation
if strcmp(method, 'ODE45')

    for ind = 1:numel(x0)
        % Create trajectory function with third order interpolation
        xTrajectory = @(t, x) uFunction(t, x, y0(ind));
        yTrajectory = @(t, y) vFunction(t, x0(ind), y);

        % Determine x and y trajectory Runge-Kutta 4th order
        % integration
        [~, x] = ode45(xTrajectory, tVec, x0(ind));
        [~, y] = ode45(yTrajectory, tVec, y0(ind));

        % assign trajectory at end time
        X(ind) = x(end);
        Y(ind) = y(end);
    end
end

% Compute trajectory using ODE45 and 3-D cubic interpolation
if strcmp(method, 'RKF45')

    [m, n] = size(x0);

    % Flatten data
    x0 = x0(:);
    y0 = y0(:);

    % RKF 45 butcher tableau

    A = [0, 2 / 9, 1 / 3, 3 / 4, 1, 5 / 6];

    B = [0, 0, 0, 0, 0, 0; ...
        2 / 9, 0, 0, 0, 0, 0; ...
        1 / 12, 1 / 4, 0, 0, 0, 0; ...
        69 / 128, -243 / 128, 135 / 64, 0, 0, 0; ...
        -17 / 12, 27 / 4, -27 / 5, 16 / 15, 0, 0; ...
        65 / 432, -5 / 16, 13 / 16, 4 / 27, 5 / 144, 0];

    CH = [47 / 450, 0, 12 / 25, 32 / 225, 1 / 30, 6 / 25];

    %     CT = [-1/150, 0, 3/100, -16/75, -1/20, 6/25];

    kx = zeros([length(x0), length(CH)]);
    ky = zeros([length(x0), length(CH)]);

    % Compute integration steps
    deltaT = gradient(tVec);

    % Loop through integration times
    for index = 1:length(tVec) - 1

        % index intergations step
        h = deltaT(index);

        % Compute integration constants
        kx(:, 1) = h * uFunction(tVec(index)*ones(size(x0))+ ...
            A(1)*h, ...
            x0, ...
            y0);
        ky(:, 1) = h * vFunction(tVec(index)*ones(size(x0))+ ...
            A(1)*h, ...
            x0, ...
            y0);

        kx(:, 2) = h * uFunction(tVec(index)*ones(size(x0))+A(2)*h, ...
            x0+kx(1, :)*B(2, :)', ...
            y0+ky(1, :)*B(2, :)');
        ky(:, 2) = h * vFunction(tVec(index)*ones(size(x0))+A(2)*h, ...
            x0+kx(1, :)*B(2, :)', ...
            y0+ky(1, :)*B(2, :)');

        kx(:, 3) = h * uFunction(tVec(index)*ones(size(x0))+A(3)*h, ...
            x0+kx(2, :)*B(3, :)', ...
            y0+ky(2, :)*B(3, :)');
        ky(:, 3) = h * vFunction(tVec(index)*ones(size(x0))+A(3)*h, ...
            x0+kx*B(3, :)', ...
            y0+ky*B(3, :)');

        kx(:, 4) = h * uFunction(tVec(index)*ones(size(x0))+A(4)*h, ...
            x0+kx(3, :)*B(4, :)', ...
            y0+ky(3, :)*B(4, :)');
        ky(:, 4) = h * vFunction(tVec(index)*ones(size(x0))+A(4)*h, ...
            x0+kx(3, :)*B(4, :)', ...
            y0+ky(3, :)*B(4, :)');

        kx(:, 5) = h * uFunction(tVec(index)*ones(size(x0))+A(5)*h, ...
            x0+kx(4, :)*B(5, :)', ...
            y0+ky(4, :)*B(5, :)');
        ky(:, 5) = h * vFunction(tVec(index)*ones(size(x0))+A(5)*h, ...
            x0+kx(4, :)*B(5, :)', ...
            y0+ky(4, :)*B(5, :)');

        kx(:, 6) = h * uFunction(tVec(index)*ones(size(x0))+A(6)*h, ...
            x0+kx(5, :)*B(6, :)', ...
            y0+ky(5, :)*B(6, :)');
        ky(:, 6) = h * vFunction(tVec(index)*ones(size(x0))+A(6)*h, ...
            x0+kx(5, :)*B(6, :)', ...
            y0+ky(5, :)*B(6, :)');

        % Compute integration
        x0 = x0 + (CH * kx')';
        y0 = y0 + (CH * ky')';

        % Compute error estimate
        %         if errorOption == true
        %         truncationErrorX(:, index) = abs(CT(1)*k1x + ...
        %             CT(2)*k2x + CT(3)*k3x...
        %             + CT(4)*k4x + CT(5)*k5x + CT(6)*k6x);
        %         truncationErrorY(:, index) = abs(CT(1)*k1y + ...
        %             CT(2)*k2y + CT(3)*k3y...
        %             + CT(4)*k4y + CT(5)*k5y + CT(6)*k6y);
        %         end
    end

    % Assign final positions to output
    X = reshape(x0, m, n);
    Y = reshape(y0, m, n);

end
end

function [U, V] = perturbation(xVec, yVec, ti, uMesh, vMesh, x0, y0, extrapolate, ugExtrap, vgExtrap)
% Function which computes trajectory of particles in a grid by
% integrating vector field of velocity values.

% Inputs
% xVec: Vector of x grid values for vector field [n x 1] (Must be in
%   ascending order)
% yVec: Vector of y grid values for vector field [m x 1] (Must be in
%   ascending order)
% tVec: Vector of time values for vector field  [p x 1]
%   (End time can be negative or positive.)
% uMesh: Matrix of x-component of velocity [m x n x p] Meshgrid style
% vMesh: Matrix of y-component of velocity [m x n x p] Meshgrid style
% x0: Matrix of intial x-position [m x n]
% y0: Matrix of intial y-position [m x n]
% extrapolate: Logical indicating whether values outside vector domain
%   should be extrapolated
% ugExtrap: Function handle for extrapolation function for u velocities.
%    Function should be in the form u = u(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest u velocity value
% vgExtrap: Function handle for extrapolation function for v velocities.
%    Function should be in the form v = v(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest v velocity value
%
% Outputs
% U: x-component of velocity [m x n]
% V: y-component of velocity [m x n]

arguments
    xVec double{mustBeVector}
    yVec double{mustBeVector}
    ti double
    uMesh double
    vMesh double
    x0 double
    y0 double
    extrapolate(1, 1) {mustBeNumericOrLogical} = false
    ugExtrap function_handle = @NOP
    vgExtrap function_handle = @NOP
end

% Initialize length
X = zeros(size(x0));
Y = zeros(size(y0));


% NaNs propogate NanNs due to cubic interpolation! Highly suggested
% that all NaNs are removed from data prior to running this function.
% Replace NaNs with mean velocity.
uMesh(isnan(uMesh)) = 0;
vMesh(isnan(vMesh)) = 0;

%Define interpolation functions
uInterp = griddedInterpolant({xVec, yVec}, ...
    permute(uMesh, [2, 1]), ...
    'linear', 'nearest'); %changed to linear to prevent spam
vInterp = griddedInterpolant({xVec, yVec}, ...
    permute(vMesh, [2, 1]), ...
    'linear', 'nearest');


    function u = uFunction(t, x, y)

        % Interpolate to velocity find values
        % (nearest value if outside domain)
        u = uInterp(x, y);

        % Create polygon of domain
        domain = [xVec(1), yVec(1); ...
            xVec(1), yVec(end); ...
            xVec(end), yVec(end); ...
            xVec(end), yVec(1); ...
            xVec(1), yVec(1)];

        % Determine values inside domain
        [in, on] = inpolygon(x, y, domain(:, 1), domain(:, 2));

        % Conditions to be outside domain
        conditions = ~in & ~on & ~isnan(x);

        % If extrapolation is true replace values of velocity
        % Use to apply taylor frozen hypothesis - ugExtrap is the DMD mode
        if extrapolate
            u(conditions) = ugExtrap(x(conditions), y(conditions), ...
                t(conditions), u(conditions));
            % Else do not advect particles outside vector field domain
        else
            u(conditions) = 0;
        end

    end

    function v = vFunction(t, x, y)
        % Interpolate to velocity find values
        % (nearest value if outside domain)
        v = vInterp(x, y);

        % Create polygon of domain
        domain = [xVec(1), yVec(1); ...
            xVec(1), yVec(end); ...
            xVec(end), yVec(end); ...
            xVec(end), yVec(1); ...
            xVec(1), yVec(1)];

        % Determine values inside domain
        [in, on] = inpolygon(x, y, domain(:, 1), domain(:, 2));

        % Conditions to be in domain
        conditions = ~in & ~on & ~isnan(x);

        % If extrapolation is true
        if extrapolate
            v(conditions) = vgExtrap(x(conditions), y(conditions), ...
                t(conditions), v(conditions));
            % Else do not advect particle
        else
            v(conditions) = 0;
        end

    end
    
    %compute velocity
    u = uFunction(ti*ones(size(x0)), x0, y0);
    v = vFunction(ti*ones(size(x0)), x0, y0);
    U = u;
    V = v;

end




% Function calculateFTLE

function out = calculateFTLE(ix, iy, flowmap_x, flowmap_y, tSpan, ...
    sxcoor, sycoor)
% Function to compute FTLE of specfic flow map location.

% Inputs
% ix: Row index of flow map
% iy: Column index of flow map
% flowmap_x: Matrix of x position of flow map [n x m]
% flowmap_y: Matrix of y position of flow map [n x m]
% tSpan: Length of integration (scalar time units (e.g. seconds))
% sxcoor: Vector of initial x positions [n x 1]
% sycoor: Vector of initial y positions [m x 1]

% Output
% out: Finite time Lyapunov exponent

% Determine size of flow map
[nx, ny] = size(flowmap_x);

% If the location is not on the edge of the flow map
if (ix - 1) * (ix - nx) < 0 && (iy - 1) * (iy - ny) < 0
    % Compute and assemble jacobian
    dPhixdX = (flowmap_x(ix+1, iy) - flowmap_x(ix-1, iy)) ...
        / (sxcoor(ix+1) - sxcoor(ix-1));
    dPhixdY = (flowmap_x(ix, iy+1) - flowmap_x(ix, iy-1)) ...
        / (sycoor(iy+1) - sycoor(iy-1));
    dPhiydX = (flowmap_y(ix+1, iy) - flowmap_y(ix-1, iy)) ...
        / (sxcoor(ix+1) - sxcoor(ix-1));
    dPhiydY = (flowmap_y(ix, iy+1) - flowmap_y(ix, iy-1)) ...
        / (sycoor(iy+1) - sycoor(iy-1));
    A = [dPhixdX, dPhixdY; dPhiydX, dPhiydY];
    % Determine max eigenvalue of streching and compute FTLE
    delta = A' * A;
    out = sqrt(max(eig(delta)));
    % If the location is on the edge of the flow map
else
    out = 0;
end


end

function NOP(varargin)
%NOP Do nothing
%
% NOP( ... )
%
% A do-nothing function for use as a placeholder when working with
%  callbacks or function handles.

% Intentionally does nothing
end