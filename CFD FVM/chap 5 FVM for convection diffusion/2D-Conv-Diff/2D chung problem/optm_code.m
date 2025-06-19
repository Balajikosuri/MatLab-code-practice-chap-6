%% 2D TDMA for Heat Transfer in a Plate
clear; clc; close all;

%% Domain and Grid Parameters
nx = 100; ny = 100;
L = 1; B = 1;
gamma = 0.000217;

dx = L / nx; dy = B / ny;
DifCof = gamma / dy;

phi_top = 0;         % Top boundary temperature
phi_bottom = 40;     % Bottom boundary temperature

%% Initialization
Phi = zeros(nx, ny);
aw = 0.5 * ones(nx, ny);
ae = -0.5 * ones(nx, ny);
an = DifCof * ones(nx, ny);
as = DifCof * ones(nx, ny);
Sp = zeros(nx, ny);
Su = zeros(nx, ny);

%% Interior Nodes Coefficients
for i = 2:nx-1
    for j = 2:ny-1
        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j);
    end
end

%% Apply Boundary Conditions
for i = 1:nx
    for j = 1:ny
        % Top boundary
        if i == 1
            an(i,j) = 0;
            Su(i,j) = 2 * DifCof * phi_top;
            Sp(i,j) = Sp(i,j) + 2 * DifCof;
            if j == 1, aw(i,j) = 0; Sp(i,j) = Sp(i,j) - 0.5; end
            if j == ny, ae(i,j) = 0; Sp(i,j) = Sp(i,j) + 0.5; end
        end

        % Bottom boundary
        if i == nx
            as(i,j) = 0;
            Su(i,j) = 2 * DifCof * phi_bottom;
            if j == 1, aw(i,j) = 0; Sp(i,j) = Sp(i,j) - 0.5; end
            if j == ny, ae(i,j) = 0; Sp(i,j) = Sp(i,j) + 0.5; end
            Sp(i,j) = Sp(i,j) + 2 * DifCof;
        end

        % Left boundary (constant heat flux = 500 kW/m²)
        if j == 1, aw(i,j) = 0; Sp(i,j) = Sp(i,j) - 0.5; end

        % Right boundary (insulated)
        if j == ny, ae(i,j) = 0; Sp(i,j) = Sp(i,j) + 0.5; end

        % Final ap with boundary sources
        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + Sp(i,j);
    end
end

%% Solve using Row-wise TDMA
[T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi);
Phi = T;

%% Add Boundary Values to Phi
T_full = frame_phi_full(Phi, phi_top, phi_bottom);

%% Plot Results
[X, Y] = meshgrid(1:(nx+2), 1:(ny+2));
contourf(X, Y, flipud(T_full), 20, 'ShowText', 'on');
colorbar;
title('2D Temperature Distribution (^oC)', 'FontSize', 16);
xlabel('x'); ylabel('y');
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

%% TDMA Function
function [T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi)
    T = Phi; tol = 1e-3; run = 1; max_res = 10;

    while max_res > tol
        T_old = T;

        for i = 1:ny
            A = zeros(1, nx); Cdash = zeros(1, nx);
            alpha = ae(i, :); beta = aw(i, :);
            D = ap(i, :);
            C = Su(i, :);

            % Add vertical neighbors
            if i > 1, C = C + as(i, :) .* T(i-1, :); end
            if i < ny, C = C + an(i, :) .* T(i+1, :); end

            % Forward Sweep
            for j = 1:nx
                if j == 1
                    A(j) = alpha(j) / D(j);
                    Cdash(j) = C(j) / D(j);
                else
                    denom = D(j) - beta(j) * A(j-1);
                    A(j) = alpha(j) / denom;
                    Cdash(j) = (C(j) + beta(j) * Cdash(j-1)) / denom;
                end
            end

            % Backward Substitution
            for j = nx:-1:1
                if j == nx
                    T(i, j) = Cdash(j);
                else
                    T(i, j) = A(j) * T(i, j+1) + Cdash(j);
                end
            end
        end

        run = run + 1;
        max_res = max(abs(T(:) - T_old(:)));
    end
end

%% Boundary Framing Function
function Phi_full = frame_phi_full(Phi_interior, phiTop, phiBottom)
    [ny, nx] = size(Phi_interior);
    Phi_full = zeros(ny + 2, nx + 2);

    Phi_full(2:end-1, 2:end-1) = Phi_interior;
    Phi_full(1, :) = phiTop;
    Phi_full(end, :) = phiBottom;
    Phi_full(2:end-1, 1) = Phi_interior(:, 1);         % Left
    Phi_full(2:end-1, end) = Phi_interior(:, end);     % Right
end
