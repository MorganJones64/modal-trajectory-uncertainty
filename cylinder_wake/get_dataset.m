%% ========================================================================
%  PART 1: Read Flow Data from Binary File
%  ========================================================================
clear all

fid = fopen('data/fixed_cylinder_atRe100', 'r');

% Read parameters: Reynolds number and Ur
params = fscanf(fid, '%f %f', [1 2]);
Re = params(1);
Ur = params(2);

fgetl(fid); % skip blank line

% Read number of timesteps and nodes
vals = fscanf(fid, '%d %d', [1 2]);
Nt = vals(1);
N_nodes = vals(2);

fgetl(fid); % skip blank line

% Preallocate arrays
times = zeros(Nt, 1);
X = zeros(Nt, N_nodes);
Y = zeros(Nt, N_nodes);
U = zeros(Nt, N_nodes);
V = zeros(Nt, N_nodes);
P = zeros(Nt, N_nodes);

% Read time series data
for n = 1:Nt
    times(n) = fscanf(fid, '%f', 1);
    
    for k = 1:N_nodes
        data = fscanf(fid, '%f %f %f %f %f', 5);
        X(n, k) = data(1);
        Y(n, k) = data(2);
        U(n, k) = data(3);
        V(n, k) = data(4);
        P(n, k) = data(5);
    end
end

fclose(fid);

% Save raw flow data
save('flow_data.mat', 'Re', 'Ur', 'times', 'X', 'Y', 'U', 'V');

%% ========================================================================
%  PART 2: Interpolate onto Regular Subgrid
%  ========================================================================
clear all
addpath('./customcolormap/')

% Load raw flow data
data = load('flow_data.mat');
times = data.times;
dt = times(2) - times(1);
Nt = length(data.times);

% Full grid resolution for interpolation
Nx = 2000;
Ny = 2000;

% Get bounds from first timestep
x0 = data.X(1, :);
y0 = data.Y(1, :);

% Create full grid
xlin = linspace(min(x0), max(x0), Nx);
ylin = linspace(min(y0), max(y0), Ny);

% Define subdomain boundaries
tol = 1;
xlin_sub = xlin(xlin >= -20 - tol & xlin <= 20 + tol);
ylin_sub = ylin(ylin >= -7 - tol  & ylin <= 7 + tol);

% Create meshgrid for subdomain
[Xsub, Ysub] = meshgrid(xlin_sub, ylin_sub);

Nx_sub = length(xlin_sub);
Ny_sub = length(ylin_sub);

% Preallocate velocity arrays
Ug_sub = zeros(Ny_sub, Nx_sub, Nt);
Vg_sub = zeros(Ny_sub, Nx_sub, Nt);

% Interpolate scattered data onto regular grid
for n = 1:Nt
    disp(['Processing timestep ', num2str(n), '/', num2str(Nt)])
    
    % Get scattered data for current timestep (column vectors)
    x_scattered = data.X(n, :)';
    y_scattered = data.Y(n, :)';
    U_scattered = data.U(n, :)';
    V_scattered = data.V(n, :)';
    
    % Create interpolants
    Fu = scatteredInterpolant(x_scattered, y_scattered, U_scattered, 'linear', 'nearest');
    Fv = scatteredInterpolant(x_scattered, y_scattered, V_scattered, 'linear', 'nearest');
    
    % Interpolate onto subgrid
    Ug_sub(:, :, n) = Fu(Xsub, Ysub);
    Vg_sub(:, :, n) = Fv(Xsub, Ysub);
end

% Prepare variables for saving
u = Ug_sub;  % Rename to u
v = Vg_sub;  % Rename to v
xMat = Xsub; % 2D meshgrid
yMat = Ysub; % 2D meshgrid
x = xlin_sub(:); % 1D vector
y = ylin_sub(:); % 1D vector
nx = Nx_sub;
ny = Ny_sub;
nt = Nt;
dx = x(2) - x(1);
dy = y(2) - y(1);

% Save processed data
fprintf('\nSaving processed data to cylindervelocity_zen.mat...\n');
save('cylindervelocity_zen.mat', ...
    'u', 'v', ...
    'xMat', 'yMat', ...
    'x', 'y', ...
    'dx', 'dt', 'dy', ...
    'nx', 'ny', 'nt', ...
    '-v7.3');
fprintf('Done!\n');
