% Clear workspace and command window
clear; clc;

% Define symbolic variables
syms A B x y

% Constants
L = 1;  % Length of the domain
H = 1;  % Height of the domain
Nx = 10; % Number of control volumes in x-direction
Ny = 10; % Number of control volumes in y-direction
lambda1 = 2 + sqrt(sym(6));
lambda2 = 2 - sqrt(sym(6));

% General solution form (adjusted to avoid unbounded growth)
phi = A*exp(x + lambda1*y) + B*exp(x + lambda2*y) + 5;

% Boundary conditions (all four corners)
eq1 = subs(phi, [x, y], [0, 0]) == 100;  % phi(0,0) = 100
eq2 = subs(phi, [x, y], [0, 1]) == 0;    % phi(0,1) = 0
eq3 = subs(phi, [x, y], [1, 0]) == 0;    % phi(1,0) = 0
eq4 = subs(phi, [x, y], [1, 1]) == 100;  % phi(1,1) = 100

% Solve for A and B
sol = solve([eq1, eq2, eq3, eq4], [A, B]);

% Convert to numeric values
A_val = double(sol.A);
B_val = double(sol.B);
lambda1_val = double(lambda1);
lambda2_val = double(lambda2);

% Display constants
disp(['A = ', num2str(A_val)]);
disp(['B = ', num2str(B_val)]);

% Define the grid for x and y
dx = L / Nx;                % Grid spacing in x-direction
dy = H / Ny;                % Grid spacing in y-direction

Px = 1:Nx;
Xi = (Px - 0.5) * dx;      % Centered grid points in x-direction
x_analytical = [0, Xi, L];  % Full x-domain with boundary

Py = 1:Ny;
Yi = (Py - 0.5) * dy;      % Centered grid points in y-direction
y_analytical = [0, Yi, H];  % Full y-domain with boundary

% Create mesh grid for plotting (now using analytical grid)
[X, Y] = meshgrid(x_analytical, y_analytical);

% Check sizes of X and Y
disp(size(X));
disp(size(Y));

% Evaluate analytical solution at grid points
PHI = A_val * exp(X + lambda1_val * Y) + B_val * exp(X + lambda2_val * Y) + 5;

% Plot the solution
figure
surf(X, Y, PHI)
xlabel('x'); ylabel('y'); zlabel('\phi(x,y)');
title('Analytical Solution of 2D Steady Convection-Diffusion-Reaction Equation');
colorbar
shading interp
