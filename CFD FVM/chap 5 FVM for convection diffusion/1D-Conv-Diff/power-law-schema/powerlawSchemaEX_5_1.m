clc; clear; close all;

%% ------------------ Physical and Boundary Conditions ------------------
U      = 2.5;      % Convection velocity
N      = 25;       % Number of control volumes
L      = 1.0;      % Domain length
Gamma  = 0.1;      % Diffusion coefficient
rho    = 1.0;      % Density
phi_A  = 1;        % Left boundary value
phi_B  = 0;        % Right boundary value

%% ------------------ Discretization Parameters ------------------
dx     = L / N;               % Control volume size
P = 1:N;
x_P = (P - 0.5) * dx;

% Initialize solution
phi = zeros(N+1,1);
phi(1)   = phi_A;
phi(end) = phi_B;

% Flux values
F = rho * U;
D = Gamma / dx;

% Power-law scheme function handle
power_law = @(Pe) max(0, (1 - 0.1 * abs(Pe)).^5);

%% ------------------ Matrix Initialization ------------------
A = zeros(N+1, N+1);
b = zeros(N+1, 1);

% Apply boundary conditions
A(1,1)     = 1;   b(1)   = phi_A;
A(end,end) = 1;   b(end) = phi_B;

%% ------------------ Assemble Internal Nodes ------------------
for i = 2:N
    Pe = F * dx / Gamma;
    aW = D * power_law(Pe) + max(F, 0);
    aE = D * power_law(Pe) + max(-F, 0);
    aP = aW + aE;
    
    A(i, i-1:i+1) = [-aW, aP, -aE];
end

%% ------------------ Solve Linear System ------------------
phi = A \ b;

%% ------------------ Analytical Solution Function ------------------
function phi = AnalyticalSolution(varargin)
    % Parse inputs using Name-Value pairs
    p = inputParser;
    addParameter(p, 'U', 1);
    addParameter(p, 'rho', 1);
    addParameter(p, 'Gamma', 0.1);
    addParameter(p, 'L', 1);
    addParameter(p, 'x', []);
    addParameter(p, 'phi_A', 1);
    addParameter(p, 'phi_B', 0);
    parse(p, varargin{:});

    % Extract parsed variables
    U      = p.Results.U;
    rho    = p.Results.rho;
    Gamma  = p.Results.Gamma;
    L      = p.Results.L;
    x      = p.Results.x;
    phi_A  = p.Results.phi_A;
    phi_B  = p.Results.phi_B;

    % Analytical solution for 1D steady-state convection-diffusion
    exponent    = (rho * U * x) / Gamma;
    denominator = exp((rho * U * L) / Gamma) - 1;
    numerator   = exp(exponent) - 1;
    phi         = phi_A + (phi_B - phi_A) * (numerator / denominator);
end

%% ------------------ Post-processing and Plotting ------------------
x_analytical = [0, x_P, L];
phi_analytical = AnalyticalSolution('U', U, 'rho', rho, 'Gamma', Gamma, ...
                                    'L', L, 'x', x_analytical, ...
                                    'phi_A', phi_A, 'phi_B', phi_B);
phi_FVM = [phi_A; phi];

figure;
hold on;
plot(x_analytical, phi_FVM, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'DisplayName', 'FVM (Power-Law Scheme)');
plot(x_analytical, phi_analytical, 'b--', 'LineWidth', 2, ...
    'DisplayName', 'Analytical Solution');

title(sprintf('Numerical vs Analytical Solution (U = %.1f, N = %d)', U, N), 'FontSize', 14);
xlabel('Distance (m)', 'FontSize', 14);
ylabel('\phi', 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;
