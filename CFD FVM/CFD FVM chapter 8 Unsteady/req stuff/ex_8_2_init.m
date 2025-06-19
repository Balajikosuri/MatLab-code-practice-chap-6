clc; close all; clear;

% Parameters
L = 0.02;              % Plate thickness (m)
N = 5;                  % Number of spatial nodes
dx = L / N;             % Grid spacing
alpha = 1e-6;          % Thermal diffusivity (m^2/s)

% Stability limit for dt
dt_limit = (dx^2) / (2 * alpha);
fprintf('Stable time step limit: %.2e s\n', dt_limit);

% Choose time steps
dt1 = 2;  % Safe time step
dt2 = dt_limit;        % Critical time step

% Simulation times
times = [40, 80, 120]; % seconds

% Geometry and discretization

dx = L / N;      % Size of control volume
P = 1:N;         
x_P = (P - 0.5) * dx;        % Cell centers
x_analytical = [0, x_P, L];  % For plotting with boundaries


% Initialize spatial and time vectors
x = linspace(0, L, N);
% x = x_analytical

% Initial condition
T0 = 200;             % °C
T_analytical = @(x, t) analytical_solution(x, t, L, alpha, 200, 100);

for dt = [dt1, dt2]
    Nt = ceil(max(times) / dt);
    T = T0 * ones(N, 1);   % Initial temperature

    % Time loop
    for step = 1:Nt
        T_new = T;

        % Internal nodes
        for i = 2:N-1
            T_new(i) = T(i) + alpha * dt / dx^2 * (T(i+1) - 2*T(i) + T(i-1));
        end

        % Boundary conditions
        T_new(1) = T_new(2);   % Insulated (Neumann BC: dT/dx = 0)
        T_new(end) = 0;        % Dirichlet BC: T = 0

        % Update
        T = T_new;

        % Compare at required times
        current_time = step * dt;
        disp(current_time)
        for tcheck = times
            if abs(current_time - tcheck) < dt/2
                % Analytical
                T_exact = T_analytical(x, tcheck);

                % Plot comparison
                figure;
                plot(x, T, 'r-o', x, T_exact, 'b-', 'LineWidth', 1.5);
                xlabel('x (m)');
                ylabel('Temperature (°C)');
                title(sprintf('T(x,t) at t = %.0f s (dt = %.1e)', tcheck, dt));
                legend('Numerical (FVM)', 'Analytical');
                grid on;
            end
        end
    end
end

% ---------- Analytical Solution Function ----------
function T = analytical_solution(x, t, L, alpha, T0, N_terms)
    T = zeros(size(x));
    for n = 1:N_terms
        lambda_n = (2*n - 1) * pi / (2*L);
        term = ((-1)^(n+1)) / (2*n - 1) * exp(-alpha * lambda_n^2 * t) .* cos(lambda_n * x);
        T = T + term;
    end
    T = T * (4 / pi) * T0;
end
