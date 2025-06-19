clc; clear;

% Parameters
L = 0.02;             % Plate thickness [m]
k = 10;               % Thermal conductivity [W/m.K]
rho_c = 1e7;          % Volumetric heat capacity [J/m^3.K]
alpha = k / rho_c;    % Thermal diffusivity [m^2/s]
Nx = 5;               % Number of nodes
dx = L / Nx;
dt = 2;               % Time step [s]
Fo = alpha * dt / dx^2;
T_init = 200;
T_left = 0;

% Spatial grid
x = linspace(0, L, Nx+1)';

% Initial temperature
T_exp = T_init * ones(Nx+1,1);
T_imp = T_exp;

% Time loop
t_end = 120;
nt = t_end / dt;

% Matrix for implicit method
% A = zeros(Nx+1);
A = zeros(Nx);
b = zeros(Nx+1,1);
ae = k/dx;
aw = k/dx;
ap_0 = rho_c*(dx/dt);

Temp_at_0_sec = ones(Nx,1)*200;
Num_Temp = Temp_at_0_sec;

for i = 1:Nx
    if i == 1 % Left node
        ap  = ae+0+ap_0;
        A(i, i)   = ap;
        A(i, i+1) = -ae;
    elseif i == Nx % Right node
        ap  = aw+0+ap_0+2*ae;
        A(i, i-1) = -aw;
        A(i, i)   =  ap;
    else % Interior nodes
        ap  = ae+aw+ap_0;
        A(i, i-1) = -aw;
        A(i, i)   = ap;
        A(i, i+1) = -ae;
    end
end
disp(A)

% % Boundary conditions
% A(1,1) = 1;              % Dirichlet at left (T = 0)
% A(end,end) = 1 + Fo;     % Neumann at right (insulated)
% A(end,end-1) = -Fo;

% Time tracking
times = [40 ];
results = struct();

for n = 1:nt
    t_current = n * dt;

    % --- EXPLICIT METHOD ---
    T_exp_new = T_exp;
    T_exp_new(1) = T_left;
    T_exp_new(end) = T_exp(end) + Fo * (T_exp(end-1) - T_exp(end)); % Neumann

    for i = 2:Nx
        T_exp_new(i) = T_exp(i) + Fo * (T_exp(i+1) - 2*T_exp(i) + T_exp(i-1));
    end
    T_exp = T_exp_new;

    % --- IMPLICIT METHOD ---
    b(1) = T_left;
    for i = 2:Nx
        b(i) = T_imp(i);
    end
    b(end) = T_imp(end);
    T_imp = A \ b;

    % --- Analytical Solution ---
    T_ana = T_left + (T_init - T_left) * erf(x / (2 * sqrt(alpha * t_current)));

    if ismember(t_current, times)
        results.(['t' num2str(t_current)]).x = x;
        results.(['t' num2str(t_current)]).explicit = T_exp;
        results.(['t' num2str(t_current)]).implicit = T_imp;
        results.(['t' num2str(t_current)]).analytical = T_ana;

        % Plot
        figure;
        plot(x, T_ana, 'k-', 'LineWidth', 2); hold on;
        plot(x, T_exp, 'ro--', 'DisplayName', 'Explicit');
        plot(x, T_imp, 'bs--', 'DisplayName', 'Implicit');
        legend('Analytical', 'Explicit FVM', 'Implicit FVM');
        xlabel('x (m)');
        ylabel('Temperature (°C)');
        title(['Temperature Profile at t = ' num2str(t_current) ' s']);
        grid on;
    end
end
