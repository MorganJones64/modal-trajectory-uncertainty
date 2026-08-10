function [dy] = currentVEC(t,yin,A,gamma,om,eps,omega,bt,at)
x = yin(1,:,:);
y = yin(2,:,:);
u = zeros(size(x)); v = u;

a = gamma * sin(om * t);
b = 1 - 2 * a;

f = a * x.^2 + b * x;
df = 2 * a * x + b;

u = 0.2-pi * A * sin(pi * f) .* cos(pi * y);
v =  pi * A * cos(pi * f) .* sin(pi * y) .* df + (sech(bt*(x - at)).^2).*(eps*sin(omega*t));

dy = [u;v];