tspan = [0 1500];
X0 = [0 0];


[tt, X] = ode23(@(tt, X) odeFunction(tt, X), tspan, X0);

X

function dXdt = odeFunction(tt, X)
M = 1060;
E = 21060;
I = 1500;

dXdt = zeros(2, 1);

dXdt(1) = X(2);
dXdt(2) = M/(E*I) * (1 + X(2)^2)^(3 / 2);
end