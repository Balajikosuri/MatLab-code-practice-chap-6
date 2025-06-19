clear;clc; close all;
% Parameters
L = 0.02;             % Plate thickness [m]
k = 10;               % Thermal conductivity [W/m.K]
rho_c = 1e7;          % Volumetric heat capacity [J/m^3.K]
Nx = 5;               % Number of nodes
dx = L / Nx;
dt = 2;               % Time step [s]

T_init = 200;
T_left = 0;

% Spatial grid
x = linspace(0, L, Nx+1);


A = zeros(Nx);
b = zeros(Nx,1);
ae = k/dx;
aw = k/dx;
ap_0 = rho_c*(dx/dt);


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
% Display the matrix A
disp('Matrix A:');
disp(A);
% Initial temperature at t = 0
previous_Time_Step_Temp = ones(Nx, 1) * T_init;  % Initial temperature (200 for each node)

% Number of time steps
t_total = [40, 80, 120];  % time points to evaluate
steps_needed = t_total / dt;

% Store results
T_all = zeros(5, length(steps_needed));

% Time loop
for step = 1:max(steps_needed)
    b = 20000 * previous_Time_Step_Temp;      % RHS = 20000 * T^n
    T_new = A \ b;        % Solve linear system
    previous_Time_Step_Temp = T_new;        % Update for next step

    % Store solution at desired time steps
    match = find(step == steps_needed);
    if ~isempty(match)
        T_all(:, match) = T_new;
        fprintf('Time = %d s\n', step * dt);
        disp(T_new);
    end
end

