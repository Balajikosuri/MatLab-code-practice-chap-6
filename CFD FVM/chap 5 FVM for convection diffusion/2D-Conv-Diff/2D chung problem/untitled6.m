matlab
% 2 D TDMA %
% Heat transfer in a 2D plate with thickness 1cm%
% Top boundary maintainsed at 100 degree Celcius.
% Constant heat flux of 500 kW/m^2 applied at west boundary
% other two boundaries are insulated. See ﻿Versteeg and Malalasekara
% Example 7.2 for problem discription.


%%%% Main Script %%%%

% 1. Define Parameters
[nx, ny, L, B, k, q] = define_parameters();

% 2. Initialize Variables
[dx, dy, aw, ae, an, as, ap, Sp, Su, T] = initialize_variables(nx, ny, L, B, k);

% 3. Apply Boundary Conditions
[aw, ae, an, as, ap, Sp, Su] = apply_boundary_conditions(nx, ny, aw, ae, an, as, ap, Sp, Su);

% 4. Solve using TDMA
[T_interior_solved, iterations, final_max_res] = solve_tdma(nx, ny, aw, ae, an, as, ap, Su, T);

% 5. Report Convergence Status
fprintf('\n--- Solution Status ---\n');
error_threshold = 0.001; % Use the same error threshold as in solve_tdma
if final_max_res < error_threshold
    fprintf('Solution converged in %i iterations.\n', iterations);
else
    fprintf('Max. iterations (%i) reached. Max residual: %.4f\n', 500, final_max_res);
end

% 6. Prepare and Display Full Temperature Matrix (for verification/debugging)
% The TDMA solves for T on the interior nodes (nx, ny).
% We need to transpose and flip for the framing function to work with the plot orientation.
T_interior_for_framing = flipud(T_interior_solved');

disp('--- Interior Temperature Matrix (for framing) ---');
disp(T_interior_for_framing);

% 7. Frame the Full Temperature Matrix including boundaries
T_full = frame_temperature_Full(T_interior_for_framing, q, dx, k);

disp('--- Full Temperature Matrix (including boundaries) ---');
disp(T_full);

% 8. Plot Temperature Distribution
plot_temperature_distribution(T_full, nx, ny);


%%%% Functions %%%%

function [nx, ny, L, B, k, q] = define_parameters()
    % Defines the problem parameters for 2D heat transfer.
    nx = 3; % number of points in x;
    ny = 4; % number of points in y;
    L = 0.3;
    B = 0.4;
    k = 1000;
    q = 500e3;
end

function [dx, dy, aw, ae, an, as, ap, Sp, Su, T] = initialize_variables(nx, ny, L, B, k)
    % Initializes the coefficient and temperature arrays.
    % Considering equal dx,dy and cell area
    dx = L / nx;
    dy = B / ny;
    A = 0.01 * 0.1; % Assuming a fixed thickness and depth for cell area
    a = (k * A) / dx;

    aw = zeros(nx, ny) + a;
    ae = zeros(nx, ny) + a;
    an = zeros(nx, ny) + a;
    as = zeros(nx, ny) + a;
    ap = zeros(nx, ny);
    Sp = zeros(nx, ny);
    Su = zeros(nx, ny);

    T = zeros(nx, ny);
end

function [aw, ae, an, as, ap, Sp, Su] = apply_boundary_conditions(nx, ny, aw, ae, an, as, ap, Sp, Su)
    % Applies boundary conditions to the coefficient arrays.

    % Interior points : ap=aw+ae+as+an
    for i = 2:nx - 1
        for j = 2:ny - 1 % Corrected inner loop limit for interior
            ap(i, j) = aw(i, j) + ae(i, j) + as(i, j) + an(i, j);
        end
    end

    % Boundaries: ap=aw+ae+as+an-Sp
    for i = 1:nx
        if i == 1 % left boundary
            for j = 1:ny
                if j == 1
                    aw(i, j) = 0;
                    as(i, j) = 0;
                    Su(i, j) = 500;
                    Sp(i, j) = 0;

                elseif j == ny
                    aw(i, j) = 0;
                    an(i, j) = 0;
                    Su(i, j) = 2500;
                    Sp(i, j) = -20;

                else
                    aw(i, j) = 0;
                    Su(i, j) = 500;
                    Sp(i, j) = 0;
                end
            end
        elseif i == nx % right boundary
            for j = 1:ny
                if j == 1
                    ae(i, j) = 0;
                    as(i, j) = 0;
                    Su(i, j) = 0;
                    Sp(i, j) = 0;

                elseif j == ny
                    an(i, j) = 0;
                    ae(i, j) = 0;
                    Su(i, j) = 2000;
                    Sp(i, j) = -20;

                else
                    Su(i, j) = 0;
                    ae(i, j) = 0;
                    Sp(i, j) = 0;
                end
            end
        else % interior columns (not left or right boundary)
            for j = 1:ny
                if j == 1 % bottom boundary in interior columns
                    as(i, j) = 0;
                    Su(i, j) = 0;
                    Sp(i, j) = 0;

                elseif j == ny % top boundary in interior columns
                    an(i, j) = 0;
                    Su(i, j) = 2000;
                    Sp(i, j) = -20;
                end
                % Note: Interior points in interior columns already have their ap set
                % This section only adjusts boundary points in interior columns
            end
        end
        % Recalculate ap for all points in the current column after boundary adjustments
        for j = 1:ny
            ap(i, j) = aw(i, j) + ae(i, j) + as(i, j) + an(i, j) - Sp(i, j);
        end
    end
end


function [T, run, max_res] = solve_tdma(nx, ny, aw, ae, an, as, ap, Su, T_initial)
    % Solves the system of equations using TDMA.
    T = T_initial; % Start with the initial temperature guess

    alpha = zeros(nx, ny);
    beta = zeros(nx, ny);
    D = zeros(nx, ny);
    C = zeros(nx, ny);

    A = zeros(nx, ny);
    Cdash = zeros(nx, ny);

    run = 1;
    error = 0.001; % set desired convergence value
    max_res = 10; % Initialize max_res high to enter the loop
    max_iter = 500;

    while max_res > error && run <= max_iter
        t = T; % Store current T for residual calculation

        for i = 1:nx
            for j = 1:ny
                alpha(i, j) = an(i, j);
                beta(i, j) = as(i, j);
                D(i, j) = ap(i, j);

                % Calculate C based on neighboring temperatures and source terms
                if i == 1 % left boundary column
                     if nx > 1 % Ensure there is a point to the east
                         C(i, j) = (ae(i, j) * T(i + 1, j)) + Su(i, j);
                     else % Handle case of a 1D problem in x if nx=1
                         C(i, j) = Su(i, j);
                     end
                elseif i == nx % right boundary column
                    if nx > 1 % Ensure there is a point to the west
                         C(i, j) = (aw(i, j) * T(i - 1, j)) + Su(i, j);
                     else % Handle case of a 1D problem in x if nx=1
                         C(i, j) = Su(i, j);
                     end
                else % interior columns
                    C(i, j) = (aw(i, j) * T(i - 1, j)) + (ae(i, j) * T(i + 1, j)) + Su(i, j);
                end
            end

            % Forward elimination (across rows, for a fixed i)
            for j = 1:ny
                if j == 1
                    A(i, j) = (alpha(i, j) / (D(i, j)));
                    Cdash(i, j) = (C(i, j) / (D(i, j)));
                else
                    A(i, j) = (alpha(i, j) / (D(i, j) - (beta(i, j) * A(i, j - 1))));
                    Cdash(i, j) = (((beta(i, j) * Cdash(i, j - 1)) + C(i, j)) / (D(i, j) - (beta(i, j) * A(i, j - 1))));
                end
            end

            % Backward substitution (across rows, for a fixed i)
            for j = ny:-1:1
                if j == ny
                    T(i, j) = Cdash(i, j);
                else
                    T(i, j) = (A(i, j) * T(i, j + 1)) + Cdash(i, j);
                end
            end
        end

        run = run + 1;
        res = abs(t - T);
        max_res = max(res(:)); % Find max residual in the entire matrix

        % fprintf('Iteration.... %i\n', run); % Uncomment to see iteration progress
    end
end


function Tfull = frame_temperature_Full(T_interior, q, dx, k)
    % Function to frame the temperature matrix with boundary conditions
    %
    % Inputs:
    % - T_interior: Interior temperature matrix (ny x nx). Expected orientation is rows=y, cols=x.
    % - q: Heat flux (W/m²)
    % - dx: Grid spacing (m)
    % - k: Thermal conductivity (W/mK)
    %
    % Output:
    % - Tfull: Full temperature matrix including boundary values (ny+2 x nx+2)

    % Get Grid Size from the Interior Matrix
    [ny_int, nx_int] = size(T_interior);

    % Compute Left Boundary Values
    % The left boundary nodes are assumed to be half a dx away from the interior nodes
    % T_left should be a column vector corresponding to the ny interior points.
    T_left = T_interior(:, 1) + (q * dx) / (2 * k);  % ny_int x 1 column vector

    % Create Full Matrix Including Boundaries (ny+2 x nx+2)
    Tfull = zeros(ny_int + 2, nx_int + 2);

    % Fill Interior Values
    Tfull(2:end-1, 2:end-1) = T_interior;

    % Apply **Top Boundary** (Fixed at 100°C)
    Tfull(1, :) = 100;

    % Apply **Left Boundary (With Heat Flux Adjustment)**
    Tfull(2:end-1, 1) = T_left;

    % Apply **Right Boundary (Insulated: Copy Preceding Column)**
    Tfull(2:end-1, end) = Tfull(2:end-1, end-1);

    % Apply **Bottom Boundary (Insulated: Copy From Row Above)**
    Tfull(end, :) = Tfull(end-1, :);

end


function plot_temperature_distribution(T_full, nx, ny)
    % Plots the temperature distribution.

    % Framing X Y Mesh for the full grid (including boundaries)
    [X, Y] = meshgrid(0:(nx + 1), 0:(ny + 1)); % Adjust meshgrid to match T_full size

    figure; % Create a new figure
    % Use surf or contourf depending on preference. surf gives a 3D surface.
    % contourf gives filled contour lines.
    % The flipud is often needed to match the plot orientation to the grid indexing.
    contourf(X, Y, flipud(T_full), 'ShowText', 'on'); % 'ShowText','on' adds labels
    colorbar; % Add color bar
    title('Temperature distribution (^o C)','FontSize',16);
    xlabel('X Grid Points');
    ylabel('Y Grid Points');
    subtitle('2D Plate - West: heatflux ; South & East: Insulated ; Top Boundary = 100 ^o C');

    % Manually add boundary labels for clarity on the plot
     % Adjust positions based on the full grid size (0 to nx+1, 0 to ny+1)
    text(nx + 1.5, (ny + 1)/2, '   Right boundary Insulated:dT_R/dx = 0', 'Rotation', 270, 'HorizontalAlignment', 'left','VerticalAlignment','middle');
    text((nx + 1)/2, -0.5, 'Bottom boundary Insulated:dT_B/dy = 0', 'HorizontalAlignment','center','VerticalAlignment','top');
    text(-0.5, (ny + 1)/2, '   steady heat flux of 500x10^3 W/m^2', 'Rotation', 90, 'HorizontalAlignment', 'left','VerticalAlignment','middle');
    text((nx + 1)/2, ny + 1.5, 'Top Boundary = 100 ^o C', 'HorizontalAlignment','center','VerticalAlignment','bottom');


    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure window
    axis equal; % Keep aspect ratio
    grid on; % Add a grid

end