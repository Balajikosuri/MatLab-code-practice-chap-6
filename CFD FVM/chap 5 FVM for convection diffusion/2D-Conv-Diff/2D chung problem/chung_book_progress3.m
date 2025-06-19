%% 2D Steady-State Heat Conduction using TDMA (Row-wise)
clear; clc; close all;

%% --- USER INPUT ---
nx = 25;                % Number of nodes in x-direction
ny = 25;                % Number of nodes in y-direction
L = 0.2;                % Length in x-direction (m)
B = 0.04;               % Height in y-direction (m)
gamma = 0.000217;       % Thermal diffusivity (m^2/s)
phi_top = 0;            % Top wall temperature (°C)
phi_left = 0;           % Left wall temperature (°C)
phi_bottom = 40;        % Bottom wall temperature (°C)

%% --- MESH SETUP ---
dx = L / nx;
dy = B / ny;
DifCof = gamma / dy;

%% --- Initialize Coefficients ---
aw = 0.5 * ones(nx, ny);
ae = -0.5 * ones(nx, ny);
an = DifCof * ones(nx, ny);
as = DifCof * ones(nx, ny);
ap = zeros(nx, ny);
Sp = zeros(nx, ny);
Su = zeros(nx, ny);
Phi = zeros(nx, ny);

%% --- Boundary Conditions and ap Setup ---
for i = 1:nx
    for j = 1:ny
        %% Top boundary (j = ny)
        if j == ny
            an(i,j) = 0;
            Su(i,j) = 2 * DifCof * phi_top;
            Sp(i,j) = Sp(i,j) + 2 * DifCof;
        end

        %% Bottom boundary (j = 1)
        if j == 1
            as(i,j) = 0;
            Su(i,j) = 2 * DifCof * phi_bottom;
            Sp(i,j) = Sp(i,j) + 2 * DifCof;
        end

        %% Left boundary (i = 1)
        if i == 1
            aw(i,j) = 0;
            Sp(i,j) = Sp(i,j) - 0.5;
            Su(i,j) = Su(i,j) + 0;  % phi_left is 0
        end

        %% Right boundary (i = nx)
        if i == nx
            ae(i,j) = 0;
            Sp(i,j) = Sp(i,j) + 0.5;
        end

        %% Final ap
        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + Sp(i,j);
    end
end

%% --- Solve using TDMA ---
[Phi, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi);

%% --- Frame Full Grid with Boundaries ---
Phi_full = frame_phi_full(Phi, phi_top, phi_bottom, phi_left);

%% --- Display and Plot ---
disp('Final Temperature Field (Phi):');
disp(Phi);

% Plot contour
[X, Y] = meshgrid(1:(nx+2), 1:(ny+2));
figure;
contourf(X, Y, flipud(Phi_full), 20, 'ShowText', 'on');
colorbar;
title('2D Temperature Distribution (^oC)');
xlabel('X (m)');
ylabel('Y (m)');
set(gca, 'YDir', 'reverse');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Plot centerline
y_vals = linspace(0, B, ny);
center_x = round(nx / 2);
figure;
plot(y_vals, flip(Phi(center_x, :)), '--o', 'LineWidth', 2);
xlabel('Y (m)');
ylabel('Temperature (^oC)');
title('Vertical Temperature Profile (Centerline)');
grid on;

%% --- TDMA Solver Function (Row-wise) ---
function [T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi)
    T = Phi;
    run = 1;
    max_res = 10;
    tol = 1e-4;
    max_iter = 10000;
    PRINT_ROWS = false; % Set to true to print row-wise matrices

    while max_res > tol && run <= max_iter
        T_old = T;

        for j = 1:ny  % Sweep each row
            alpha = zeros(1, nx);
            beta  = zeros(1, nx);
            diag  = zeros(1, nx);
            rhs   = zeros(1, nx);
            A     = zeros(1, nx);
            Cdash = zeros(1, nx);

            for i = 1:nx
                alpha(i) = ae(i, j);
                beta(i)  = aw(i, j);
                diag(i)  = ap(i, j);
                rhs(i)   = Su(i, j);

                % Include vertical neighbors
                if j > 1
                    rhs(i) = rhs(i) + as(i, j) * T(i, j-1);
                end
                if j < ny
                    rhs(i) = rhs(i) + an(i, j) * T(i, j+1);
                end
            end

            % Print system matrix (optional)
            if PRINT_ROWS
                fprintf('\nRow #%d (TDMA):\n', j);
                fprintf('aw:   '); fprintf('%8.4f ', beta); fprintf('\n');
                fprintf('ap:   '); fprintf('%8.4f ', diag); fprintf('\n');
                fprintf('ae:   '); fprintf('%8.4f ', alpha); fprintf('\n');
                fprintf('RHS:  '); fprintf('%8.4f ', rhs); fprintf('\n');
            end

            % TDMA forward sweep
            for i = 1:nx
                if i == 1
                    A(i) = alpha(i) / diag(i);
                    Cdash(i) = rhs(i) / diag(i);
                else
                    denom = diag(i) - beta(i) * A(i-1);
                    if abs(denom) < 1e-12
                        error('Zero denominator in TDMA at i=%d, j=%d', i, j);
                    end
                    A(i) = alpha(i) / denom;
                    Cdash(i) = (rhs(i) + beta(i) * Cdash(i-1)) / denom;
                end
            end

            % TDMA backward substitution
            for i = nx:-1:1
                if i == nx
                    T(i, j) = Cdash(i);
                else
                    T(i, j) = A(i) * T(i+1, j) + Cdash(i);
                end
            end
        end

        run = run + 1;
        max_res = max(abs(T(:) - T_old(:)));
    end
end

%% --- Add Boundary Rows/Columns to Phi ---
function Phi_full = frame_phi_full(Phi_interior, phiTop, phiBottom, phiLeft)
    [nx, ny] = size(Phi_interior);
    Phi_full = zeros(nx+2, ny+2);
    Phi_full(2:end-1, 2:end-1) = Phi_interior;
    Phi_full(1, :) = phiTop;
    Phi_full(end, :) = phiBottom;
    Phi_full(2:end-1, 1) = phiLeft;
    Phi_full(2:end-1, end) = Phi_interior(:, end);
end
