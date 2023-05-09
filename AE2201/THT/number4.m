%% Main Program
clc, clear

% Input the NIM
nim = char(input("Insert your NIM: ", "s"));

% Parse NIM for information
global E M I;
E = str2double(nim(4:8));
M = str2double(nim(5:8));
I = 1500;

% Get boundary conditions
x0 = 0;
y0 = [0, 0];

% Get other information
L = 200;
h = 20;

[x2, y2] = modifiedEulerRK2nd(@number4ODE, x0, y0, L, h);
[x4, y4] = classicalRK4th(@number4ODE, x0, y0, L, h);

fig2 = figure();
plot(x2, y2(:,2), "Marker","+")
hold on
grid on
plot(x4, y4(:,2), "Marker","o")
xlim([0 L])
ylim([min(y4(:,1)) max(y4(:,2))])
legend(["Modified Euler Runge-Kutta 2nd Order" "Classical Runge-Kutta 4th Order"])
xlabel("x (mm)")
ylabel("z")
title(sprintf("Beam Deflection with Bernoulli-Euler Theory\nM = %d kg mm, E = %d kg/mm^{2}, I = %d mm^{4}, h = %d mm", M, E, I, h))
hold off

fig1 = figure();
plot(x2, y2(:,1), "Marker","+")
hold on
grid on
plot(x4, y4(:,1), "Marker","o")
xlim([0 L])
ylim([min(y4(:,1)) max(y4(:,1))])
legend(["Modified Euler Runge-Kutta 2nd Order" "Classical Runge-Kutta 4th Order"])
xlabel("x (mm)")
ylabel("y (mm)")
title(sprintf("Beam Deflection with Bernoulli-Euler Theory\nM = %d kg mm, E = %d kg/mm^{2}, I = %d mm^{4}, h = %d mm", M, E, I, h))
hold off

%% Modified Euler Runge-Kutta 2nd Order Function
function [xSol, ySol] = modifiedEulerRK2nd(diffEq, x0, y0, L, h)
n = floor(L / h);
order = length(y0);

xSol = zeros(n + 1, 1);
ySol = zeros(n + 1, order);

xSol(1) = x0;
ySol(1,:) = y0;

i = 2;
x = x0;
y = y0;

while x < L
    h = min(h, L - x);
    K1 = diffEq(x, y);
    K2 = diffEq(x + h, y + K1 * h);
    x = x + h;
    y = y + (K1 + K2) * h/2;
    xSol(i) = x;
    ySol(i,:) = y;
    i = i + 1;
end
end

%% Classical Runge-Kutta 4th Order
function [xSol, ySol] = classicalRK4th(diffEq, x0, y0, L, h)
n = floor(L / h);
order = length(y0);

xSol = zeros(n + 1, 1);
ySol = zeros(n + 1, order);

xSol(1) = x0;
ySol(1,:) = y0;

i = 1;
x = x0;
y = y0;

while x < L
    i = i + 1;
    h = min(h, L - x);
    K1 = diffEq(x, y);
    K2 = diffEq(x + h / 2, y + K1 * h / 2);
    K3 = diffEq(x + h / 2, y + K2 * h / 2);
    K4 = diffEq(x + h, y + K3 * h);
    x = x + h;
    y = y + (K1 + 2*K2 + 2*K3 + K4) * h/6;
    xSol(i) = x;
    ySol(i,:) = y;
end
end

%% Differential Function
function yz__1 = number4ODE(x, y)
global M E I;

z = y(1,2);

y__1 = z;
z__1 = M/(E*I) * (1 + z^2)^(3/2);

yz__1 = [y__1, z__1];
end
