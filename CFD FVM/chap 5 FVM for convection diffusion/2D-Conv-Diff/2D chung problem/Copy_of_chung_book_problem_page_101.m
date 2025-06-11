%% 2D TDMA for Steady-State Heat Conduction
clear; clc; close all;

%% Problem Setup
nx = 3;             % Number of points in x-direction
ny = 3;              % Number of points in y-direction

ny = ny + 1 ;           %We are calculating the Discritization equation for this problem
                      % right wa as aP 
L = 0.2;              % Length in x-direction
B = 0.04;             % Height in y-direction
gamma = 0.000217;     % Thermal diffusivity
phi_top = 0;          % Top boundary temperature
phi_left = 0;         % left boundary temperature
phi_bottom = 40;      % Bottom boundary temperature

dx = L / nx;
dy = B / ny;
DifCof = gamma / dy;

%% Initialize Coefficients
aw = ones(nx, ny) * 0.5;
ae = ones(nx, ny) * -0.5;
an = ones(nx, ny) * DifCof;
as = ones(nx, ny) * DifCof;
ap = zeros(nx, ny);
Sp = zeros(nx, ny);
Su = zeros(nx, ny);
Phi = zeros(nx, ny);

%% Compute ap for interior points
for i = 2:nx-1
    for j = 2:ny-1
        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j);
    end
end

%% Boundary Conditions 
for i = 1:nx
    for j = 1:ny
        % Top boundary
        if i == 1
            % fprintf('balaji   (i,j) = (%f,%f)\n',i,j);
            an(i,j) = 0;
            Su(i,j) = 2 * DifCof * phi_top;
            if j == 1
                aw(i,j) = 0;
                Sp(i,j) = 2 * DifCof - 0.5;
                % Su(i,j) = 2 * DifCof * phi_top;
                    
            elseif j==ny
                ae(i,j) = 0;
                Sp(i,j) = 2 * DifCof + 0.5;
                % Su(i,j) = 2 * DifCof * phi_top;
            else
                Sp(i,j) = 2 * DifCof;
            end
        end
       
        %% Bottom boundary
        if i == nx
            as(i,j) = 0;
            Su(i,j) = 2 * DifCof * phi_bottom;
            if  j == 1
                aw(i,j) = 0;
                Sp(i,j) = 2 * DifCof - 1/2;
                % fprintf('balaji line 39 (i,j) = (%f,%f)\n',i,j);
            elseif j == ny
                ae(i,j) = 0;
                Sp(i,j) = 2 * DifCof + 1/2;
            else
                Sp(i,j) = 2 * DifCof;
            end
        end
        %% Left boundary
        if j == 1
            % fprintf('balaji line 39 (i,j) = (%f,%f)\n',i,j);
            aw(i,j) = 0;
            Sp(i,j) = Sp(i,j) - 0.5;
            if i == 1 
                an(i,j) = 0;
                Sp(i,j) = 2 * DifCof - 0.5;
                Su(i,j) = 2 * DifCof * phi_top;
            elseif i == nx
                % fprintf('85 balaji   (i,j) = (%f,%f)\n',i,j);
                as(i,j) = 0;
                Sp(i,j) = 2 * DifCof - 0.5;
                Su(i,j) = 2 * DifCof * phi_bottom;
            else
                % fprintf('balaji   (i,j) = (%f,%f)\n',i,j);
                Sp(i,j) = -1/2;
                Su(i,j) = 0;
            end
            
        end
        % Right boundary
        if j == ny
            if i == 1

            elseif i == nx
                %d
            else
                fprintf('balaji (i,j) = (%f,%f)\n',i,j);
            end

            ae(i,j) = 0;
            Sp(i,j) = 0.5;
        end

        % Final ap after adding Sp
        ap(i,j) = aw(i,j) + ae(i,j) + as(i,j) + an(i,j) + Sp(i,j);
    end
end

%% Solve using TDMA (row-wise in x-direction)
[T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi);
Phi = T;
% disp(Phi);

%% Frame with boundaries
phiTop = phi_top;
phiBottom = phi_bottom;
T_full = frame_phi_full(Phi, phiTop, phiBottom);

%% Plot Contour
[X, Y] = meshgrid(1:(nx+2), 1:(ny+2));
contourf(X, Y, flipud(T_full'), 'ShowText', 'on'); 
set(gca, 'YDir', 'reverse');               % ✅ flip y-axis for correct orientation
colorbar;
title('Temperature Distribution (^oC)', 'FontSize', 22);
xlabel('X (grid points)');
ylabel('Y (grid points)');

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

%% Central Line Plot (vertical profile at center x)
y = linspace(0, B, ny);
center_x = round(nx / 2);
figure;
plot(y, flip(Phi(center_x, :)), '--o', 'LineWidth', 2);
xlabel('y (m)');
ylabel('Temperature (^oC)');
title('Temperature along Centerline in Y-direction');
grid on;

%% --- TDMA Solver Function ---
function [T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi)
    T = Phi;
    run = 1;
    max_res = 10;
    tol = 1e-3;

    while max_res > tol
        T_old = T;

        for i = 1:ny  % For each row (fixed y)
            alpha = zeros(1, nx);
            beta = zeros(1, nx);
            D = zeros(1, nx);
            C = zeros(1, nx);
            A = zeros(1, nx);
            Cdash = zeros(1, nx);

            for j = 1:nx
                alpha(j) = ae(j, i);
                beta(j) = aw(j, i);
                D(j) = ap(j, i);
                C(j) = Su(j, i);

                % Add contributions from vertical neighbors
                if i > 1
                    C(j) = C(j) + as(j, i) * T(j, i - 1);
                end
                if i < ny
                    C(j) = C(j) + an(j, i) * T(j, i + 1);
                end
            end

            % TDMA forward sweep
            for j = 1:nx
                if j == 1
                    A(j) = alpha(j) / D(j);
                    Cdash(j) = C(j) / D(j);
                else
                    denom = D(j) - beta(j) * A(j - 1);
                    A(j) = alpha(j) / denom;
                    Cdash(j) = (C(j) + beta(j) * Cdash(j - 1)) / denom;
                end
            end

            % TDMA backward substitution
            for j = nx:-1:1
                if j == nx
                    T(j, i) = Cdash(j);
                else
                    T(j, i) = A(j) * T(j + 1, i) + Cdash(j);
                end
            end
        end

        run = run + 1;
        max_res = max(abs(T(:) - T_old(:)));
    end
end

%% --- Function to Frame the Interior Matrix with Boundaries ---
function Phi_full = frame_phi_full(Phi_interior, phiTop, phiBottom)
    [ny, nx] = size(Phi_interior);
    Phi_full = zeros(ny + 2, nx + 2);
    Phi_full(2:end-1, 2:end-1) = Phi_interior;
    Phi_full(1, :) = phiTop;
    Phi_full(end, :) = phiBottom;
    Phi_full(2:end-1, 1) = Phi_interior(:, 1);
    Phi_full(2:end-1, end) = Phi_interior(:, end);
end
