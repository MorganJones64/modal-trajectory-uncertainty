clear all
fid = fopen('fixed_cylinder_atRe100','r');

%%

% First line: Re, Ur
params = fscanf(fid, '%f %f', [1 2]);
Re = params(1);
Ur = params(2);

fgetl(fid); % skip blank line

% Nt and N_nodes
vals = fscanf(fid, '%d %d', [1 2]);
Nt = vals(1);
N_nodes = vals(2);

fgetl(fid); % skip blank line

times = zeros(Nt,1);
X = zeros(Nt, N_nodes);
Y = zeros(Nt, N_nodes);
U = zeros(Nt, N_nodes);
V = zeros(Nt, N_nodes);
P = zeros(Nt, N_nodes);

for n = 1:Nt
    times(n) = fscanf(fid, '%f', 1);

    for k = 1:N_nodes
        data = fscanf(fid, '%f %f %f %f %f', 5);
        X(n,k) = data(1);
        Y(n,k) = data(2);
        U(n,k) = data(3);
        V(n,k) = data(4);
        P(n,k) = data(5);
    end
end

fclose(fid);

save('flow_data.mat',"Re","Ur","times","X","Y","U","V")
%%
clear all
addpath('./customcolormap/')

data = load('flow_data.mat');
times = data.times
dt = times(2)-times(1)
Nt = length(data.times);

% Grid resolution
Nx = 2000;
Ny = 2000;

% --- Build FULL grid ONCE (only for bounds) ---
x0 = data.X(1,:);
y0 = data.Y(1,:);

xlin = linspace(min(x0), max(x0), Nx);
ylin = linspace(min(y0), max(y0), Ny);
tol = 1
% --- Define subdomain directly ---
xlin_sub = xlin(xlin >= -20 -tol & xlin <= 20 + tol);
ylin_sub = ylin(ylin >= -7 - tol  & ylin <= 7 + tol);

[Xsub, Ysub] = meshgrid(xlin_sub, ylin_sub);

Nx_sub = length(xlin_sub);
Ny_sub = length(ylin_sub);

% --- Preallocate storage ---
Ug_sub = zeros(Ny_sub, Nx_sub, Nt);
Vg_sub = zeros(Ny_sub, Nx_sub, Nt);

% --- Time loop ---
for n = 1:Nt
    disp(['Processing timestep ', num2str(n), '/', num2str(Nt)])

    x = data.X(n,:)';   % column vectors
    y = data.Y(n,:)';
    U = data.U(n,:)';
    V = data.V(n,:)';

    % Interpolants
    Fu = scatteredInterpolant(x, y, U, 'linear', 'nearest');
    Fv = scatteredInterpolant(x, y, V, 'linear', 'nearest');

    % ✅ Interpolate DIRECTLY onto subgrid (better)
    Ug_s = Fu(Xsub, Ysub);
    Vg_s = Fv(Xsub, Ysub);

    % Store
    Ug_sub(:,:,n) = Ug_s;
    Vg_sub(:,:,n) = Vg_s;
end

% --- Save everything ---
save('cylindervelocity_zen.mat', ...
    'Ug_sub', 'Vg_sub', ...
    'xMat', 'yMat', ...
    'x', 'y', ...
    'data','dx','dt','dy','nx','ny','nt','', '-v7.3')

%x should be a vec y should be a vec, xMat is Xsub, yMat is Ysub
\
