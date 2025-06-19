clc; clear; close all

%% Domain Setup
L = 1.0;                % x-length (streamwise)
h = 0.04;               % y-height
Nx = 501;
Ny = 41;

dx = L / (Nx - 1);
dy = h / (Ny - 1);
x = linspace(0, L, Nx);
y = linspace(0, h, Ny);

%% Physical Properties
u = 1.0;                % Velocity in x-direction (m/s)
v = 0.0;                % Velocity in y-direction
nu = 0.000217;          % Kinematic viscosity

%% --- Péclet Number Calculation ---
Pe_x = u * L / nu;
Pe_y = (v * h) / nu;

fprintf('Péclet Number (x-direction): Pe_x = %.2f → %s\n', Pe_x, dominanceCategory(Pe_x));
fprintf('Péclet Number (y-direction): Pe_y = %.2f → %s\n\n', Pe_y, dominanceCategory(Pe_y));

%% Initialize Field (e.g. φ could be temperature or velocity)
phi = zeros(Ny, Nx);
phi(1,1) = 40;          % Bottom inlet energized
phi(2:end,1) = 0;       % Rest of inlet

%% Marching in x-direction (explicit, like time-stepping)
for i = 1:Nx-1
    for j = 2:Ny-1
        % Central difference for vertical diffusion (pure)
        phi_yy = (phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) / dy^2;

        % March in x (explicit scheme)
        phi(j,i+1) = phi(j,i) + dx/u * nu * phi_yy;
    end

    % Reapply boundary conditions
    phi(1,i+1) = 40;
    phi(end,i+1) = 0;
end

%% Plot Results at Specific x Locations
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

%% --- Péclet Categorization Function ---
function category = dominanceCategory(Pe)
    if Pe < 0.1
        category = 'Diffusion-dominated';
    elseif Pe > 10
        category = 'Convection-dominated';
    else
        category = 'Mixed regime';
    end
end
