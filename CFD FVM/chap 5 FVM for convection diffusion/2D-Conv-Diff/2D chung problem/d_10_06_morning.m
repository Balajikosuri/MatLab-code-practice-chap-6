clc; clear; close all

% Domain
L = 1.0;
h = 0.04;
Nx = 501;%501
Ny = 41;%41

dx = L / (Nx - 1);
dy = h / (Ny - 1);
x = linspace(0, L, Nx);
y = linspace(0, h, Ny);

% Physical properties
u = 1;                 
nu = 0.000217;

% Initialize field
phi = zeros(Ny, Nx);

% Boundary conditions
% phi(1,:) = 0;      % Top wall
phi(end,:) = 40;     % Bottom wall
phi(2:end,1) = 0;   % left wall inlet (all zero)
% phi(1,1) = 8;      % Inlet bottom edge energized
disp(phi)

% Marching in x-direction
for i = 1:Nx-1
    for j = 2:Ny-1
        % Diffusion in y-direction
        phi_yy = (phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) / dy^2;
        phi(j,i+1) = phi(j,i) + dx/u * nu * phi_yy;
    end

    % Reapply wall BCs
    phi(1,i+1) = 40;
    phi(end,i+1) = 0;
end

% Plot results at selected x locations
x_sample = [0.2, 0.4, 0.6, 1.0];
colors = 'rbgk';
figure; hold on
for k = 1:length(x_sample)
    idx = round(x_sample(k) / dx) + 1;
    plot(phi(:, idx), y, colors(k), 'DisplayName', ['x = ' num2str(x_sample(k))])
end
xlabel('\phi'); ylabel('y'); grid on
title('Convection-Diffusion: \phi(y) at Selected x')
legend show
