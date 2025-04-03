clc; clear; close all;
%  2D diffusion equation using the line-by-line TDMA (Thomas Algorithm).

%% Grid and Parameters
Lx = 1;  % Length in x-direction
Ly = 1;  % Length in y-direction
Nx = 20; % Number of grid points in x
Ny = 20; % Number of grid points in y
dx = Lx/Nx;
dy = Ly/Ny;
tol = 1e-5; % Convergence tolerance
max_iter = 5000; % Maximum iterations

% Diffusion coefficient (constant)
gama = 1; 

% Coefficients for the discretized equation
    % Along X axis
aW = gama/dx;
aE = gama/dx;

    % Along Y axis
aS = gama/dy;
aN = gama/dy;
aP = -2 * (aW + aS); % Central coefficient

%% Initialize Field Variable (phi)
phi = zeros(Nx, Ny); % Initial guess

% Boundary Conditions
phi(:, 1) = 400;  % Left boundary (Dirichlet)
phi(:, end) =300; % Right boundary (Dirichlet)
phi(1, :) = 200;  % Bottom boundary
phi(end, :) = 500; % Top boundary

%% Iterative TDMA Solution
for iter = 1:max_iter
    phi_old = phi;
    
    % Sweep in x-direction (solving for y-lines)
    for j = 2:Ny-1
        % Setup tridiagonal system
        A = diag(aP * ones(Nx-2,1)) + diag(aE * ones(Nx-3,1), 1) + diag(aW * ones(Nx-3,1), -1);
        d = -aS * phi(2:Nx-1, j-1) - aN * phi(2:Nx-1, j+1);
        d(1) = d(1) - aW * phi(1, j); % Left boundary
        d(end) = d(end) - aE * phi(Nx, j); % Right boundary
        
        % Solve using Thomas Algorithm (TDMA)
        phi(2:Nx-1, j) = tdma_solver(A, d);
    end

    % Sweep in y-direction (solving for x-lines)
    for i = 2:Nx-1
        % Setup tridiagonal system
        A = diag(aP * ones(Ny-2,1)) + diag(aN * ones(Ny-3,1), 1) + diag(aS * ones(Ny-3,1), -1);
        d = -aW * phi(i-1, 2:Ny-1)' - aE * phi(i+1, 2:Ny-1)';
        d(1) = d(1) - aS * phi(i, 1); % Bottom boundary
        d(end) = d(end) - aN * phi(i, Ny); % Top boundary
        
        % Solve using Thomas Algorithm (TDMA)
        phi(i, 2:Ny-1) = tdma_solver(A, d)';
    end

    % Check Convergence
    error = max(max(abs(phi - phi_old)));
    if error < tol
        fprintf('Converged in %d iterations\n', iter);
        break;
    end
end

%% Plot Results
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
figure;
contourf(X, Y, phi', 20, 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('2D Diffusion Solution using Line-by-Line TDMA');

%% TDMA Solver Function
function x = tdma_solver(A, d)
    % Solves Ax = d using the Thomas Algorithm (TDMA)
    n = length(d);
    a = diag(A, -1);
    b = diag(A);
    c = diag(A, 1);
    
    % Forward Elimination
    for i = 2:n
        w = a(i-1) / b(i-1);
        b(i) = b(i) - w * c(i-1);
        d(i) = d(i) - w * d(i-1);
    end
    
    % Back Substitution
    x = zeros(n,1);
    x(n) = d(n) / b(n);
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end
end
