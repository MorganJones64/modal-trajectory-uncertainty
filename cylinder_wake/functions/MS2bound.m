function [cseIntegral, deltaInfty, xPos, yPos] = MS2bound(ufMesh, vfMesh, ugMesh, vgMesh, xVec, yVec, ...
    t0, tLength, tStep, dt, filename, xMinROI, xMaxROI, ...
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
% uExtrap: Function handle for extrapolation function for u velocities.
%    Function should be in the form u = u(x, y, t, nearest) where x, y,
%    and t are the time and location of the particle outside the domain and
%    nearest is the closest u velocity value
% vExtrap: Function handle for extrapolation function for v velocities.
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
    ufMesh single
    vfMesh single
    ugMesh single
    vgMesh single
    xVec double{mustBeVector}
    yVec double{mustBeVector}
    t0(1, 1) double
    tLength(1, 1) double
    tStep(1, 1) double
    dt(1, 1) double
    filename(1,:) char
    xMinROI(1, 1) double = min(xVec)
    xMaxROI(1, 1) double = max(xVec)
    yMinROI(1, 1) double = max(yVec)
    yMaxROI(1, 1) double = max(yVec)
    nx(1, 1) double{mustBeInteger} = 100
    ny(1, 1) double{mustBeInteger} = 100
    method(1, :) char ...
        {mustBeMember(method, {'Euler', 'RK4'})} = 'RK4'
    options.xMask double{mustBeVector}
    options.yMask double{mustBeVector}
    options.extrap(1, 1) {mustBeNumericOrLogical} = false
    options.ufExtrap function_handle
    options.vfExtrap function_handle
    options.ugExtrap function_handle
    options.vgExtrap function_handle
end

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

    % Define time vector
    %tVec = (0:sign(tStep):tStep) * dt; % Issue with extrapolation 0-dt is fixed
    dt = sign(tStep)*dt;
    tend = tLength + t0 - sign(tStep);
    tLoop = t0:tStep:tend;
    % Initalize position matrix
    xPos = repmat(xPos, 1, 1, length(tLoop));
    yPos = repmat(yPos, 1, 1, length(tLoop));
    
    cseIntegral = zeros(length(tLoop),nx, ny);
    deltaInfty = zeros(length(tLoop),nx, ny); 
    sigma_cse = zeros(abs(tLength),nx, ny);
    %loop through end periods
    for t = t0:tStep:tend
        t
        for s = t0:tStep:t %sweeping parameter (s)
            %s
                sVec1 = t0:tStep:s;
                sVec2 = (s+sign(tStep)):tStep:t;
            %xs = xPos(50,50,2,1)
            for i = sVec1
                %t
                % Determine the index
                tIndex = find(sVec1 == i);
                if sVec1 == t0
                    continue
                end
                % Compute trajectories FIX tVec
                 [xPos(:, :, tIndex+1), yPos(:, :, tIndex+1)] = ...
                    trajectory(xVec, yVec, i, dt, ...
                    ufMesh(:, :, i), ...
                    vfMesh(:, :, i), ...
                    xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
                    options.extrap, @options.ufExtrap, @options.vfExtrap);
            end
    
            % save g value at time s sVec1(end)
            [gx0(:,:,abs(s-t0)+1), gy0(:,:,abs(s-t0)+1)] = ...
                 perturbation(xVec,yVec,i,...
                 ugMesh(:, :, i),...
                 vgMesh(:, :, i),...
                 xPos(:, :, tIndex), yPos(:, :, tIndex),...
                 options.extrap, @options.ugExtrap, @options.vgExtrap);
            
            %Advect from s->t
            for i = sVec2
                if sVec2 == t0
                    continue
                end
                %t
                % Determine the index
                tIndex = find(sVec2 == i);
                % Compute trajectories FIX tVec
                [xPos(:, :, tIndex+1), yPos(:, :, tIndex+1)] = ...
                    trajectory(xVec, yVec, i, dt, ...
                    ufMesh(:, :, i), ...
                    vfMesh(:, :, i), ...
                    xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
                    options.extrap, @options.ufExtrap, @options.vfExtrap);
            end
            %xe = xPos(50,50,2,end)
            % Determine components of zJacobian. szcoord will only have 3 elements,
            % based on three plane locations
            [dPhiXdY, dPhiXdX] = gradient(xPos(:, :, end), sycoor, sxcoor);
    
            [dPhiYdY, dPhiYdX] = gradient(yPos(:, :, end), sycoor, sxcoor);
        
            % NEEDS WORK
            % Construct Jacobian for each element
            A1 = cat(3, dPhiXdX, dPhiXdY);
    
            A2 = cat(3, dPhiYdX, dPhiYdY);
        
            % Create 5-D matrix of jacobians
            % nx x ny x nz x (3 x 3) -> dPhid
            B = cat(4, A1, A2);
        
            % Change order of dimensions for page transpose
            B = permute(B, [3, 4, 1, 2]);

            % Compute stretching for only the middle z plane (index 2)
            delta = pagemtimes(B, pagetranspose(B));

            %max(B,[],'all')
            % Compute FTLE for each point
            for j = 1:numel(sigma_cse(1,:))
                sigma_cse(abs(s-t0)+1,j) = sqrt(max(eig(delta(:, :, j))));
            end
        end
            % If extrapolation is false stop integrating trajectories when
            % particle exits vector field domain  
        %for the specific end interval of integration t, compute the CSE across the whole grid
        % *** NEEDS TO BE FIXED FOR ALL OTHER CASES
        for j = 1:numel(sigma_cse(1,:))
            cseIntegral(abs(t-t0)+1,j) = trapz(sigma_cse(:,j)*abs(dt));
        end
        % *** NEEDS TO BE FIXED FOR ALL OTHER CASES
        deltaInfty(abs(t-t0)+1,:,:) = max(sqrt(gx0.^2 + gy0.^2),[],3);
        %ss=t0-fstart+1
        
        clear gx0 gy0
        if mod(t,20) == 0 
            save(filename,'cseIntegral','deltaInfty','t','-append')
        end
    end
    %reformat to nx x ny x nt
    deltaInfty = permute(deltaInfty,[2,3,1]);
    cseIntegral = permute(cseIntegral,[2,3,1]);

else
    % % Find bounds of vector field data
    % xMin = min(xVec);
    % xMax = max(xVec);
    % yMin = min(yVec);
    % yMax = max(yVec);
    % 
    % % Define the initial region of the fluid to track
    % sxcoor = linspace(xMinROI, xMaxROI, nx);
    % sycoor = linspace(yMinROI, yMaxROI, ny);
    % 
    % % Initial position
    % [xPos, yPos] = ndgrid(sxcoor, sycoor);
    % 
    % % Check for mask
    % if isfield(options, 'xMask') && isfield(options, 'yMask')
    %     % find positions inside mask and set to nan
    %     [in, ~] = inpolygon(xPos, yPos, options.xMask, options.yMask);
    %     xPos(in) = nan;
    %     yPos(in) = nan;
    % end
    % 
    % 
    % % Define time vector
    % %tVec = (0:sign(tStep):tStep) * dt; % Issue with extrapolation 0-dt is fixed
    % dt = sign(tStep)*dt;
    % tend = tLength + t0 - sign(tStep);
    % tLoop = t0:tStep:tend;
    % % Initalize position matrix
    % xPos = repmat(xPos, 1, 1, length(tLoop));
    % yPos = repmat(yPos, 1, 1, length(tLoop));
    % 
    % %gives [nt,nx,ny]
    % xPos = repmat(reshape(xPos, [1, nx, ny]), length(tLoop), 1, 1);
    % yPos = repmat(reshape(yPos, [1, nx, ny]), length(tLoop), 1, 1);
    % 
    % cseIntegral = zeros(length(tLoop),nx, ny);
    % deltaInfty = zeros(length(tLoop),nx, ny); 
    % sigma_cse = zeros(abs(tLength),nx, ny);
    % 
    % % Initialize flags
    % calcsigma = ones(nx, ny, 'logical');
    % calcdelta = ones(nx, ny, 'logical');
    % outDomain = zeros(nx, ny, 'logical');
    % 
    % %loop through end periods
    % for t = t0:tStep:tend
    %     t
    %     for s = t0:tStep:t %sweeping parameter (s)
    %         %s
    %             sVec1 = t0:tStep:s;
    %             sVec2 = (s+sign(tStep)):tStep:t;
    %         %xs = xPos(50,50,2,1)
    %         for i = sVec1
    %             %t
    %             % Determine the index
    %             tIndex = find(sVec1 == i);
    %             if sVec1 == t0
    %                 continue
    %             end
    %             % Compute trajectories FIX tVec
    %              [xPos(tIndex+1, ~outDomain), yPos(tIndex+1, ~outDomain)] = ...
    %                 trajectory(xVec, yVec, i, dt, ...
    %                 ufMesh(:, :, i), ...
    %                 vfMesh(:, :, i), ...
    %                 xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
    %                 options.extrap, @options.ufExtrap, @options.vfExtrap);
    % 
    %         end
    % 
    %         % save g value at time s sVec1(end)
    %         [gx0(:,:,abs(s-t0)+1), gy0(:,:,abs(s-t0)+1)] = ...
    %              perturbation(xVec,yVec,i,...
    %              ugMesh(:, :, i),...
    %              vgMesh(:, :, i),...
    %              xPos(:, :, tIndex), yPos(:, :, tIndex),...
    %              options.extrap, @options.ugExtrap, @options.vgExtrap);
    % 
    %         %Advect from s->t
    %         for i = sVec2
    %             if sVec2 == t0
    %                 continue
    %             end
    %             %t
    %             % Determine the index
    %             tIndex = find(sVec2 == i);
    %             % Compute trajectories FIX tVec
    %             [xPos(:, :, tIndex+1), yPos(:, :, tIndex+1)] = ...
    %                 trajectory(xVec, yVec, i, dt, ...
    %                 ufMesh(:, :, i), ...
    %                 vfMesh(:, :, i), ...
    %                 xPos(:, :, tIndex), yPos(:, :, tIndex), method, ...
    %                 options.extrap, @options.ufExtrap, @options.vfExtrap);
    %         end
    %         %xe = xPos(50,50,2,end)
    %         % Determine components of zJacobian. szcoord will only have 3 elements,
    %         % based on three plane locations
    %         [dPhiXdY, dPhiXdX] = gradient(xPos(:, :, end), sycoor, sxcoor);
    % 
    %         [dPhiYdY, dPhiYdX] = gradient(yPos(:, :, end), sycoor, sxcoor);
    % 
    %         % NEEDS WORK
    %         % Construct Jacobian for each element
    %         A1 = cat(3, dPhiXdX, dPhiXdY);
    % 
    %         A2 = cat(3, dPhiYdX, dPhiYdY);
    % 
    %         % Create 5-D matrix of jacobians
    %         % nx x ny x nz x (3 x 3) -> dPhid
    %         B = cat(4, A1, A2);
    % 
    %         % Change order of dimensions for page transpose
    %         B = permute(B, [3, 4, 1, 2]);
    % 
    %         % Compute stretching for only the middle z plane (index 2)
    %         delta = pagemtimes(B, pagetranspose(B));
    % 
    %         %max(B,[],'all')
    %         % Compute FTLE for each point
    %         for j = 1:numel(sigma_cse(1,:))
    %             sigma_cse(abs(s-t0)+1,j) = sqrt(max(eig(delta(:, :, j))));
    %         end
    %     end
    %         % If extrapolation is false stop integrating trajectories when
    %         % particle exits vector field domain  
    %     %for the specific end interval of integration t, compute the CSE across the whole grid
    %     % *** NEEDS TO BE FIXED FOR ALL OTHER CASES
    %     for j = 1:numel(sigma_cse(1,:))
    %         cseIntegral(abs(t-t0)+1,j) = trapz(sigma_cse(:,j)*abs(dt));
    %     end
    %     % *** NEEDS TO BE FIXED FOR ALL OTHER CASES
    %     deltaInfty(abs(t-t0)+1,:,:) = max(sqrt(gx0.^2 + gy0.^2),[],3);
    %     %ss=t0-fstart+1
    % 
    %     clear gx0 gy0
    %     if mod(t,20) == 0 
    %         save(filename,'cseIntegral','deltaInfty','t','-append')
    %     end
    % end
    % %reformat to nx x ny x nt
    % deltaInfty = permute(deltaInfty,[2,3,1]);
    % cseIntegral = permute(cseIntegral,[2,3,1]);
        
        % % Find bounds of vector field data
        %  xMin = min(xVec);
        % xMax = max(xVec);
        % yMin = min(yVec);
        % yMax = max(yVec);
        % 
        % % Define the initial region of the fluid to track
        % sxcoor = linspace(xMinROI, xMaxROI, nx);
        % sycoor = linspace(yMinROI, yMaxROI, ny);
        % 
        % % Initial position
        % [xPosTemp, yPosTemp] = ndgrid(sxcoor, sycoor);
        % 
        % % Check for mask
        % if isfield(options, 'xMask') && isfield(options, 'yMask')
        %     % find positions inside mask and set to nan
        %     [in, ~] = inpolygon(xPosTemp, yPosTemp, ...
        %         options.xMask, options.yMask);
        %     xPosTemp(in) = nan;
        %     yPosTemp(in) = nan;
        % else
        %     % If no mask exists set all points to be outside mask
        %     in = zeros(size(xPosTemp));
        % end
    
        % % Initialize FTLE value
        % sigma = zeros(nx, ny);
        % 
        % % Integration time length
        % tSpan = abs(tLength) * dt;
        % 
        % % Initialize flags
        % calcFTLE = ones(nx, ny, 'logical');
        % outDomain = zeros(nx, ny, 'logical');
        % 
        % % Create time vector for loop
        % tLoop = t0:tStep:t0 + tLength;
        % 
        % % Initalize position matrix
        % xPos = repmat(zeros(nx, ny), 1, 1, length(tLoop));
        % yPos = repmat(zeros(nx, ny), 1, 1, length(tLoop));
        % 
        % for t = tLoop
        % 
        %     tIndex = find(tLoop == t);
        % 
        %     xPos(:, :, tIndex) = xPosTemp;
        %     yPos(:, :, tIndex) = yPosTemp;
        % 
        %     % If the point is inside, calculate it's trajectory
        %     [xPosTemp(~outDomain), yPosTemp(~outDomain)] = ...
        %         trajectory(xVec, yVec, tVec, ...
        %         ufMesh(:, :, t:sign(tStep)*1:t+tStep), ...
        %         vfMesh(:, :, t:sign(tStep)*1:t+tStep), ...
        %         xPosTemp(~outDomain), yPosTemp(~outDomain), method);
        % 
        % 
        %     % If the mask exists
        %     if isfield(options, 'xMask') && isfield(options, 'yMask')
        %         % Check for positions inside mask
        %         [in, ~] = inpolygon(xPosTemp, yPosTemp, ...
        %             options.xMask, options.yMask);
        %     end
        % 
        %     % Determine particles that have trajectories outside the domain
        %     index = (xPosTemp - xMin) .* (xPosTemp - xMax) >= 0 | ...
        %         (yPosTemp - yMin) .* (yPosTemp - yMax) >= 0 | in;
        % 
        %     % Change flag to indicate they have departed
        %     outDomain(index) = true;
        % 
        %     % Find adjacent points
        %     M = zeros(size(xPosTemp));
        %     M(index) = 1;
        %     index = conv2(M, [1, 1, 1; 1, 1, 1; 1, 1, 1], 'same') > 0;
        % 
        %     % Find row and column of points
        %     [ix, iy] = find(index);
        % 
        %     % Loop through points that have exited domain and their neighbors
        %     for i = 1:numel(ix)
        %         % Compute FTLE
        %         sigma(ix(i), iy(i)) = calculateFTLE(ix(i), iy(i), ...
        %             xPosTemp, yPosTemp, ...
        %             tSpan, sxcoor, sycoor);
        %     end
        % 
        %     % Change flag to indicate these points FTLEs have been computed
        %     calcFTLE(index) = false;
        % end
        % 
        % % calculate FTLE for all the points inside the domain at the last
        % % time step.
        % 
        % % Find row and column of points
        % [ix, iy] = find(calcFTLE);
        % 
        % % Loop through points that have exited domain and their neighbors
        % for i = 1:numel(ix)
        %     % Compute FTLE
        %     sigma(ix(i), iy(i)) = calculateFTLE(ix(i), iy(i), ...
        %         xPosTemp, yPosTemp, ...
        %         tSpan, sxcoor, sycoor);
        % end
end
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

function NOP(varargin)
%NOP Do nothing
%
% NOP( ... )
%
% A do-nothing function for use as a placeholder when working with
%  callbacks or function handles.

% Intentionally does nothing
end