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
T_left = 0;           % Left boundary condition (not applied directly in code)

% Spatial grid
x = linspace(0, L, Nx+1);  % x from 0 to L (6 points)
x_inner = x(2:end);        % Interior nodes (excluding boundary at x=0)

% Matrix A and RHS b
A = zeros(Nx);
b = zeros(Nx, 1);
ae = k / dx;
aw = k / dx;
ap_0 = rho_c * (dx / dt);

for i = 1:Nx
    if i == 1 % Left node
        ap = ae + ap_0;
        A(i, i)   = ap;
        A(i, i+1) = -ae;
    elseif i == Nx % Right node
        ap = aw + ap_0 + 2*ae;
        A(i, i-1) = -aw;
        A(i, i)   = ap;
    else % Interior nodes
        ap = ae + aw + ap_0;
        A(i, i-1) = -aw;
        A(i, i)   = ap;
        A(i, i+1) = -ae;
    end
end

% Initial condition
previous_T = ones(Nx, 1) * T_init;

% Time points to evaluate
t_total = [40, 80, 120];
steps_needed = t_total / dt;
T_all = zeros(Nx, length(steps_needed));

% Time loop
for step = 1:max(steps_needed)
    b = 20000 * previous_T;
    T_new = A \ b;
    previous_T = T_new;

    % Store results at selected times
    match = find(step == steps_needed);
    if ~isempty(match)
        T_all(:, match) = T_new;
    end
end

%% Analytical Solution Function
analytical_T = zeros(Nx+1, length(t_total));
n_terms = 20000;  % Number of terms in Fourier series

for idx = 1:length(t_total)
    t = t_total(idx);
    for xi = 1:length(x)
        sum_series = 0;
        for n = 1:n_terms
            lambda_n = (2*n - 1)*pi / (2*L);
            coeff = ((-1)^(n+1)) / (2*n - 1);
            sum_series = sum_series + coeff * exp(-alpha * lambda_n^2 * t) * cos(lambda_n * x(xi));
        end
        analytical_T(xi, idx) = (4*T_init/pi) * sum_series;
    end
end
% add the boundary nodes for 
left_Bc = T_all(2, :);
right_Bc = zeros(1, size(T_all, 2));

T_Full = [left_Bc;T_all,right_Bc];
disp(T_Full);

disp(analytical_T);

P = 1:Nx;
x_P = (P - 0.5) * dx;
x_analytical = [0, x_P, L];

%% Plotting
figure;
for idx = 1:length(t_total)
    subplot(1, 3, idx);
    plot(x, analytical_T(:, idx), 'o-', 'DisplayName', 'Analytical');
    hold on;
    plot(x_analytical, T_all(:, idx), 'x--', 'DisplayName', 'Numerical');
    xlabel('x [m]');
    ylabel('Temperature [°C]');
    title(['Time = ', num2str(t_total(idx)), ' s']);
    legend;
    grid on;
end
%%