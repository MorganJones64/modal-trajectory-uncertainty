 %Duffing model
% dx/dt = y
% dy/dt = gamma*cos(omega*t)-delta*y-alpha*x-beta*x^3;

function dx = gVEC(t,yin,omega,bt,at)
x = yin(1,:,:);
u = zeros(size(x));
y = yin(2,:,:);
v = zeros(size(y));
% x is a two dimensional state vector
% x(1),x(2),x(3) - x, y
% dx(1),dx(2),dx(3) - dx/dt, dy/dt, dz/dt
u=u;
v=v + (sech(bt*(x-at)).^2).*(sin(omega*t));
dx = [u;v];