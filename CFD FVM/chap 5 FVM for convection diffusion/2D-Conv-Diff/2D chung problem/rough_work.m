clear; close all; clc;

%% Grid and Parameters
nx = 500; % Points in x
ny = 40;  % Points in y
L = 0.2;
B = 0.04;

gamma = 0.000217;
phi_top = 0;
phi_bottom = 40;

dx = L / nx;
dy = B / ny;
DifCof = gamma / dy;

%% Coefficient Arrays — now sized [ny, nx]
aw = ones(ny, nx) * (1/2);
ae = ones(ny, nx) * (-1/2);
an = ones(ny, nx) * DifCof;
as = ones(ny, nx) * DifCof;
ap = zeros(ny, nx);
Sp = zeros(ny, nx);
Su = zeros(ny, nx);
Phi = zeros(ny, nx);

%% Apply interior and boundary coefficients
for j = 1:ny
    for i = 1:nx
        % Boundaries
        if i == 1     % Left
            aw(j,i) = 0;
            Sp(j,i) = Sp(j,i) - 1/2;
        elseif i == nx  % Right
            ae(j,i) = 0;
            Sp(j,i) = Sp(j,i) + 1/2;
        end
        if j == 1     % Top
            an(j,i) = 0;
            Su(j,i) = Su(j,i) + 2 * DifCof * phi_top;
            Sp(j,i) = Sp(j,i) + 2 * DifCof;
        elseif j == ny % Bottom
            as(j,i) = 0;
            Su(j,i) = Su(j,i) + 2 * DifCof * phi_bottom;
            Sp(j,i) = Sp(j,i) + 2 * DifCof;
        end
        % Total ap
        ap(j,i) = aw(j,i) + ae(j,i) + an(j,i) + as(j,i) + Sp(j,i);
    end
end

% Avoid division by zero
ap(ap == 0) = 1e-12;

%% TDMA Solver (row-by-row in y)
[T, run, max_res] = solve_tdma_rowwise(ny, nx, aw, ae, an, as, ap, Su, Phi);
Phi = T;

%% Pad with boundaries for full plot
T_full = frame_phi_full(Phi, phi_top, phi_bottom);
[X, Y] = meshgrid(1:(nx+2), 1:(ny+2));

%% Plot
figure;
contourf(X, Y, flipud(T_full), 20, 'ShowText', 'on');
colorbar;
title('Temperature distribution (^oC)', 'FontSize', 20);
xlabel('X'); ylabel('Y');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

%% =================== FUNCTIONS ===================

function [T, run, max_res] = solve_tdma_rowwise(ny, nx, aw, ae, an, as, ap, Su, Phi)
    T = Phi;
    run = 1;
    max_res = 10;
    tol = 1e-3;

    while max_res > tol
        T_old = T;
        for j = 1:ny
            alpha = ae(j, :);
            beta = aw(j, :);
            D = ap(j, :);
            C = Su(j, :);
            if j > 1
                C = C + as(j, :) .* T(j - 1, :);
            end
            if j < ny
                C = C + an(j, :) .* T(j + 1, :);
            end

            A = zeros(1, nx);
            Cdash = zeros(1, nx);

            for i = 1:nx
                if i == 1
                    denom = D(i);
                else
                    denom = D(i) - beta(i) * A(i - 1);
                end
                if abs(denom) < 1e-12
                    denom = 1e-12;
                end
                A(i) = alpha(i) / denom;
                if i == 1
                    Cdash(i) = C(i) / denom;
                else
                    Cdash(i) = (C(i) + beta(i) * Cdash(i - 1)) / denom;
                end
            end

            for i = nx:-1:1
                if i == nx
                    T(j,i) = Cdash(i);
                else
                    T(j,i) = A(i) * T(j,i + 1) + Cdash(i);
                end
            end
        end
        run = run + 1;
        max_res = max(abs(T(:) - T_old(:)));
    end
end

function Phi_full = frame_phi_full(Phi_interior, phiTop, phiBottom)
    [ny, nx] = size(Phi_interior);
    Phi_full = zeros(ny + 2, nx + 2);
    Phi_full(2:end-1, 2:end-1) = Phi_interior;
    Phi_full(1, :) = phiTop;
    Phi_full(end, :) = phiBottom;
    Phi_full(2:end-1, 1) = Phi_interior(:, 1);
    Phi_full(2:end-1, end) = Phi_interior(:, end);
end
