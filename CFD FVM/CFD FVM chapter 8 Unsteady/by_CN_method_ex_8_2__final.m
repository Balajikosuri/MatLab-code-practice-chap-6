clear; clc; close all;

% Parameters.
Nx = 5;             % Number of control volumes
dt = 2;               % Time step [s]
L = 0.02;             % Plate thickness [m]
k = 10;               % Thermal conductivity [W/m.K]
rho_c = 1e7;          % Volumetric heat capacity [J/m^3.K]
alpha = k / rho_c;    % Thermal diffusivity
dx = L / Nx;
%The dt Limit by the Fully Implicit Scheme
% Compute stability limit for Explicit FVM
dt_lim = (rho_c * dx^2) / (k);

fprintf('The dt Limit = %.4f by the Fully Implicit Scheme\n', dt_lim);

if dt >= dt_lim
    % Popup warning and stop execution
    msg = sprintf('ERROR:\nYour dt = %.4f exceeds the stability limit (%.4f) for the Fully Implicit method.', dt, dt_lim);
    errordlg(msg, 'Time Step Too Large');
    return;
else
    disp('The FVM Solution is Given By:');
end




T_init = 200;         % Initial temperature
T_Right = 0;          % Dirichlet boundary at right

% Spatial grid
P = 1:Nx;
x_P = (P - 0.5) * dx;       % control volume centers
x_analytical = [0, x_P, L]; % for full analytical comparison

% Coefficients
ae = k / dx;
aw = k / dx;
ap_0 = rho_c * dx / dt;

% Matrix A
A = zeros(Nx);
for i = 1:Nx
    if i == 1
        ap = ae/2 + ap_0;
        A(i, i) = ap;
        A(i, i+1) = -ae/2;
    elseif i == Nx
        ap = aw/2 + ap_0 + 2*ae/2; % Neumann-like correction
        A(i, i-1) = -aw/2;
        A(i, i) = ap;
    else
        ap = ae/2 + aw/2 + ap_0;
        A(i, i-1) = -aw/2;
        A(i, i) = ap;
        A(i, i+1) = -ae/2;
    end
end

% Initial condition
T = ones(Nx, 1) * T_init;

% Time loop and plotting
t_total = [40, 80, 120];
% t_total = [40];
steps_needed = floor(t_total / dt);
% steps_needed = (t_total / dt);

T_storage = zeros(Nx+2, length(t_total));  % include boundary nodes

figure;
% for each_time = t_total

% for step = 1:max(steps_needed)
for step = 1:max(steps_needed)
    b = ones(Nx,1);
    for i = 1:Nx
        if i == 1
            ap_0_coeff = (ap_0-ae/2);
            b(i) = ap_0_coeff*T(i) + (ae/2)*T(i+1);
        elseif i == Nx
            ap_0_coeff = (ap_0-aw/2-ae);
            b(i) = (aw/2)*T(i-1) + ap_0_coeff*T(i)+ae*T_Right;
        else
            ap_0_coeff = (ap_0-ae/2-aw/2);
            b(i) = (aw/2)*T(i-1) + ap_0_coeff*T(i) + (ae/2)*T(i+1);
        end

    end
    T = A \ b;

    % Plot and store at selected times
    t_current = step * dt;
    idx_plot = find(t_current == t_total);
    if ~isempty(idx_plot)
        % Analytical solution
        T_analytical = zeros(size(x_analytical));
        for xi = 1:length(x_analytical)
            sum_series = 0;
            for n = 1:20000
                lambda_n = (2*n - 1)*pi / (2*L);
                coeff = ((-1)^(n+1)) / (2*n - 1);
                sum_series = sum_series + coeff * exp(-alpha * lambda_n^2 * t_current) * cos(lambda_n * x_analytical(xi));
            end
            T_analytical(xi) = (4*T_init/pi) * sum_series;
        end

        % Store full temperature: [boundary-left estimate; T; boundary-right]
        T_full = [T(1); T; T_Right];
        T_storage(:, idx_plot) = T_full;

        % Plot
        subplot(1, 3, idx_plot);
        plot(x_analytical, T_analytical, 'o-', 'LineWidth', 1.5);
        hold on;
        plot(x_analytical, T_full, 'rx--', 'LineWidth', 1.5);
        xlabel('x [m]');
        ylabel('Temperature [°C]');
        title(['t = ', num2str(t_current), ' s']);
        legend('Analytical', 'Numerical');
        grid on;
    end
end
% end
%% Print Comparison Table 
fprintf('\n%-6s | %-7s | %-17s | %-17s | %-10s\n', ...
        'Time', 'Node', 'Numerical (°C)', 'Analytical (°C)', '% Error');
fprintf(repmat('-', 1, 65)); fprintf('\n');

for t_idx = 1:length(t_total)
    t_val = t_total(t_idx);
    T_num = T_storage(:, t_idx);  % including boundaries

    % Analytical recomputation at x_analytical
    T_ana_vals = zeros(length(x_analytical), 1);
    for xi = 1:length(x_analytical)
        sum_series = 0;
        for n = 1:20000
            lambda_n = (2*n - 1)*pi / (2*L);
            coeff = ((-1)^(n+1)) / (2*n - 1);
            sum_series = sum_series + coeff * exp(-alpha * lambda_n^2 * t_val) * cos(lambda_n * x_analytical(xi));
        end
        T_ana_vals(xi) = (4*T_init/pi) * sum_series;
    end

    for node = 1:length(x_analytical)
        num_val = T_num(node);
        ana_val = T_ana_vals(node);
        err = ((num_val - ana_val) / ana_val) * 100;
        fprintf('%-6d | %-7d | %-17.2f | %-17.2f | %-10.2f\n', ...
                t_val, node-1, num_val, ana_val, err);
    end
    fprintf(repmat('-', 1, 65)); fprintf('\n');
end

