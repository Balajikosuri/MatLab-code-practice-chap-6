clc; clear; close all

%% Grid and Physical Parameters
L = 1.0;
h = 0.04;
Nx = 436;           % Cells in x suitable Nx are form 436 
Ny = 41;            % Cells in y use 40 due to Inital phi 40;

dx = L / (Nx - 1);  % Cell size in x
dy = h / (Ny - 1);  % Cell size in y
x = linspace(0, L, Nx);
y = linspace(0, h, Ny);

u = 1.0;            % Flow speed in x (m/s)
nu = 0.000217;      % Kinematic viscosity

%% FVM Initialization
phi = zeros(Ny, Nx);    % φ(j, i) — i: x-direction, j: y-direction

% Apply Boundary Conditions
phi(1,1) = 40;          % Bottom-left corner (inlet bottom)
phi(2:end,1) = 0;       % All other inlet values

%% Marching in x using FVM logic
for i = 1:Nx-1
    for j = 2:Ny-1
        % Finite volume central difference in y (diffusion)
        dphi_dy2 = (phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) / dy^2;

        % Update using explicit forward step in x
        phi(j,i+1) = phi(j,i) + (dx/u) * nu * dphi_dy2;
    end

    % Reapply BCs at each x-slice
    phi(1,i+1) = 40;   % Bottom wall
    phi(end,i+1) = 0;  % Top wall
end

%% Plot Comparison (with FTCS targets)
x_sample = [0.2, 0.4, 0.6, 1.0];
colors = 'rbgk';
figure; hold on
for k = 1:length(x_sample)
    idx = round(x_sample(k)/dx) + 1;
    plot(phi(:, idx), y, colors(k), 'DisplayName', ['x = ' num2str(x_sample(k))])
end
xlabel('\phi'); ylabel('y'); grid on
title('FVM: \phi(y) at selected x-locations (Verification)')
legend show

%% Optional: Print numerical diffusion number
d = nu * dx / (u * dy^2);
fprintf('Numerical diffusion number (d) = %.4f\n', d);
