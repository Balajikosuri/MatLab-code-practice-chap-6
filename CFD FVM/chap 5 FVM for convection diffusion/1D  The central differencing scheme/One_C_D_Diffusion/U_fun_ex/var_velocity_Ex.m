clc; clear; close all;

% --- Parameters ---
L = 1.0;         % Length of domain (m)
Gamma = 0.1;     % Diffusion coefficient (kg/m/s)
rho = 1.0;       % Density (kg/m³)
u0 = 1.0;        % Base velocity (m/s)
alpha = 1.0;     % Linear growth factor for u(x)
N = 5;          % Number of control volumes
dx = L / N;      % Width of control volume

% Discretize space
x_nodes = ((1:N) - 0.5) * dx;   % Cell center positions
phi_A = 1;                      % Left BC
phi_B = 0;                      % Right BC

% Initialize matrix system
A = zeros(N, N);
B = zeros(N, 1);

% Loop over all nodes
for i = 1:N
    % Cell center
    x_i = x_nodes(i);

    % Face positions
    x_W = x_i - dx / 2;
    x_E = x_i + dx / 2;

    % Velocity at faces
    u_W = u0 * (1 + alpha * x_W / L);
    u_E = u0 * (1 + alpha * x_E / L);
    u_P = u0 * (1 + alpha * x_i / L);

    % Mass flux at faces
    F_W = rho * u_W;
    F_E = rho * u_E;
    F = rho*u_P;

    % Diffusion coefficients
    D_W = Gamma / dx;
    D_E = Gamma / dx;

    if i == 1
        % Left boundary node
        a_E = D_E - F_E / 2;
        a_P = a_E + -2*D_W-F;
        A(i, i) = a_P;
        A(i, i+1) = -a_E;
        B(i) = (2*D_W+F) * phi_A;

    elseif i == N
        % Right boundary node
        a_W = D_W + F_W / 2;
        a_P = a_W + -2*D_E +F;
        A(i, i) = a_P;
        A(i, i-1) = -a_W;
        B(i) = (2*D_E - F) * phi_B;

    else
        % Internal node
        a_W = D_W + F_W / 2;
        a_E = D_E - F_E / 2;
        a_P = a_W + a_E;
        A(i, i) = a_P;
        A(i, i-1) = -a_W;
        A(i, i+1) = -a_E;
        B(i) = 0;
    end
end

% Solve the system
phi_numerical = A \ B;

% --- Analytical Solution (Approximate using Integration) ---
% Define u(x), and solve ODE: d/dx(rho u(x) phi) = d/dx(Gamma dphi/dx)
% For variable u(x), use integrating factor method

% Create fine mesh for analytical solution
P = 1:N; % Node indices
x_P = (P - 0.5) * dx;  % Apply x_P = x_0 + (P - 1/2) * dx

% Include boundary points
x_analytical = [0, x_P, L];
u_x = u0 * (1 + alpha * x_analytical / L);
Pe_x = rho * u_x / Gamma;

% Compute integrating factor:
% φ(x) = 1 - ∫₀ˣ μ(s) ds / ∫₀ᴸ μ(s) ds
mu = @(s) exp(rho / Gamma * (u0 * (s + alpha * s.^2 / (2*L))));
integral_mu = cumtrapz(x_analytical, mu(x_analytical));
total_mu = integral_mu(end);
phi_analytical = 1 - integral_mu / total_mu;

phi_FVM_Num_full = [phi_A;phi_numerical;phi_B];
% --- Plot ---
figure;
plot(x_analytical, phi_FVM_Num_full, 'ro-', 'LineWidth', 2, 'MarkerSize', 8); hold on;
plot(x_analytical, phi_analytical, 'b-', 'LineWidth', 2);
xlabel('\xi = (P - 0.5)\Deltax', 'FontSize', 14);
ylabel('\phi(x)', 'FontSize', 14);
title(['FVM vs Analytical | u(x) = u_0(1 + \alpha x/L),  \alpha = ', num2str(alpha)], 'FontSize', 14);
legend('Numerical (FVM)', 'Analytical');
grid on;
