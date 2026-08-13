function [dy] = crossflowVEC(t,yin,U,A,gamma,omega_gyre,eps,omega_cross,bt,at)
x = yin(1,:,:);
y = yin(2,:,:);
u = zeros(size(x)); v = u;

a = gamma * sin(omega_gyre * t);
b = 1 - 2 * a;

f = a * x.^2 + b * x;
df = 2 * a * x + b;

u = U - pi * A * sin(pi * f) .* cos(pi * y);
v =  pi * A * cos(pi * f) .* sin(pi * y) .* df + (sech(bt*(x - at)).^2).*(eps*sin(omega_cross * t));

dy = [u;v];