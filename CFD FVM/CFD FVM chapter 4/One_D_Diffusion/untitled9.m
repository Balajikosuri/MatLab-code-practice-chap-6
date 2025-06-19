clear; clc; close all;

% Domain and grid
nx = 3;  % Including boundaries
ny = 4;
L = 0.2;
B = 0.04;
gamma = 0.000217;

phi_top = 0;
phi_left = 0;
phi_bottom = 40;
phi_right = 0;  % Dirichlet or Neumann (adjust below)

dx = L / (nx - 1);
dy = B / (ny - 1);
DifCof = gamma / dy;

% Initialize fields
Phi = zeros(nx, ny);  % [i,j] = [x, y]

% Coefficients
aw = ones(nx, ny) * 0.5;
ae = ones(nx, ny) * -0.5;
an = ones(nx, ny) * DifCof;
as = ones(nx, ny) * DifCof;
ap = zeros(nx, ny);
Su = zeros(nx, ny);
Sp = zeros(nx, ny);

% === Apply Boundary Conditions ===

% Top boundary (j = 1)
Phi(:, 1) = phi_top;

% Bottom boundary (j = ny)
Phi(:, ny) = phi_bottom;

% Left boundary (i = 1)
Phi(1, :) = phi_left;

% Right boundary (i = nx)
Phi(nx, :) = phi_right;  % OR: Phi(nx, :) = Phi(nx-1, :) for insulated

% === Compute Coefficients and Source Terms ===

for i = 2:nx-1
    for j = 2:ny-1
        ap(i,j) = aw(i,j) + ae(i,j) + an(i,j) + as(i,j);
    end
end

% === Gauss-Seidel Solver ===
tol = 1e-4;
max_iter = 10000;
iter = 0;
max_res = 10;

while max_res > tol
    Phi_old = Phi;
    

    for i = 2:nx-1
        for j = 2:ny-1
            Phi(i,j) = ( ...
                aw(i,j) * Phi(i-1,j) + ...
                ae(i,j) * Phi(i+1,j) + ...
                as(i,j) * Phi(i,j-1) + ...
                an(i,j) * Phi(i,j+1) + ...
                Su(i,j)) / ap(i,j);
        end
    end

    % Optional: re-enforce BCs (in case of accidental overwrite)
    Phi(:, 1) = phi_top;
    Phi(:, ny) = phi_bottom;
    Phi(1, :) = phi_left;
    Phi(nx, :) = phi_right;

    max_res = max(abs(Phi(:) - Phi_old(:)));
end

fprintf('Converged in %d iterations. Max residual = %.5e\n', iter, max_res);
disp('Final temperature matrix:')
disp(Phi);

% === Plot ===
[X, Y] = meshgrid(0:dx:L, 0:dy:B);
figure;
contourf(X, Y, Phi', 20, 'ShowText', 'on');
set(gca, 'YDir', 'normal');
colorbar;
title('Steady-State Temperature (Gauss-Seidel)');
xlabel('x (m)');
ylabel('y (m)');
