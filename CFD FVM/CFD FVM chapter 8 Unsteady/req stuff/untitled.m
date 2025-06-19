clear; clc; close all;

% Parameters
L = 0.02;             % Plate thickness [m]
k = 10;               % Thermal conductivity [W/m.K]
rho_c = 1e7;          % Volumetric heat capacity [J/m^3.K]
alpha = k / rho_c;    % Thermal diffusivity
Nx = 5;               % Number of nodes
dx = L / Nx;
dt = 2;               % Time step [s]
T_init = 200;         % Initial temperature
T_right = 0;           % Dirichlet boundary (east side)

% Spatial grid
x = linspace(0, L, Nx+1);   % 6 points (including boundary)
x_inner = x(2:end);         % Interior nodes (excluding 0 m)

% Coefficients
ae = k / dx;
aw = k / dx;
ap_0 = rho_c * dx / dt;

% Matrix A
A = zeros(Nx);
for i = 1:Nx
    if i == 1
        ap = ae + ap_0;
        A(i, i) = ap;
        A(i, i+1) = -ae;
    elseif i == Nx
        ap = aw + ap_0 + 2*ae; % Adjusted for Neumann
        A(i, i-1) = -aw;
        A(i, i) = ap;
    else
        ap = ae + aw + ap_0;
        A(i, i-1) = -aw;
        A(i, i) = ap;
        A(i, i+1) = -ae;
    end
end

% Initial condition
T = ones(Nx, 1) * T_init;

% Time loop and plotting
t_total = [40, 80, 120];
steps_needed = t_total / dt;
x_plot = [x(2), x(3), x(4), x(5), x(6)];

figure;
for step = 1:max(steps_needed)
    b = 20000 * T;
    T = A \ b;

    % Plot at selected times
    t_current = step * dt;
    idx_plot = find(t_current == t_total);
    if ~isempty(idx_plot)
        % Analytical solution
        T_analytical = zeros(size(x));
        for xi = 1:length(x)
            sum_series = 0;
            for n = 1:20000
                lambda_n = (2*n - 1)*pi / (2*L);
                coeff = ((-1)^(n+1)) / (2*n - 1);
                sum_series = sum_series + coeff * exp(-alpha * lambda_n^2 * t_current) * cos(lambda_n * x(xi));
            end
            T_analytical(xi) = (4*T_init/pi) * sum_series;
        end
        T_full = [T(1);T;T_right];
        subplot(1, 3, idx_plot);
        plot(x, T_analytical, 'o-', 'LineWidth', 1.5); hold on;
        plot(x_plot, T_full, 'rx--', 'LineWidth', 1.5);
        xlabel('x [m]');
        ylabel('Temperature [°C]');
        title(['t = ', num2str(t_current), ' s']);
        legend('Analytical', 'Numerical');
        grid on;
    end
end

