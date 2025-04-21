
clc; clear; close all;

% Physical and boundary values
U = 2.5;
N = 25;
L = 1.0;
Gamma = 0.1;
rho = 1.0;
phi_A = 1;
phi_B = 0;


% Discretization
dx = L/N;
% Initialize solution array
phi = zeros(N+1,1);
phi(1) = phi_A;
phi(end) = phi_B;

% Flux coefficients
F = rho * U;
D = Gamma / dx;

% Power-law scheme function
power_law = @(Pe) max(0, (1 - 0.1*abs(Pe)).^5);

% Matrix initialization
A = zeros(N+1,N+1);
b = zeros(N+1,1);

% Apply boundary conditions
A(1,1) = 1;
b(1) = phi_A;
A(end,end) = 1;
b(end) = phi_B;

% Assemble equations for internal nodes
for i = 2:N
    Pe = F*dx/Gamma;
    aW = D*power_law(Pe) + max(F,0);
    aE = D*power_law(Pe) + max(-F,0);
    aP = aW + aE;
    
    A(i,i-1:i+1) = [-aW, aP, -aE];
end

% Solve linear system
phi = A\b;
function phi = AnalyticalSolution(varargin)
    % Parse Name-Value pairs
    p = inputParser;

    addParameter(p, 'U', 1);               % Velocity
    addParameter(p, 'rho', 1);            % Density
    addParameter(p, 'Gamma', 0.1);        % Diffusion coefficient
    addParameter(p, 'L', 1);              % Domain length
    addParameter(p, 'x', []);             % x locations where phi is evaluated
    addParameter(p, 'phi_A', 1);          % Left boundary value
    addParameter(p, 'phi_B', 0);          % Right boundary value

    parse(p, varargin{:});
    U      = p.Results.U;
    rho    = p.Results.rho;
    Gamma  = p.Results.Gamma;
    L      = p.Results.L;
    x      = p.Results.x;
    phi_A  = p.Results.phi_A;
    phi_B  = p.Results.phi_B;

    % Analytical formula for 1D steady convection-diffusion problem
    exponent = (rho * U * x) / Gamma;
    denominator = exp((rho * U * L) / Gamma) - 1;
    numerator = exp(exponent) - 1;

    phi = phi_A + (phi_B - phi_A) * (numerator / denominator);
end




% Discretization
P = 1:N;
x_P = (P - 0.5) * dx;
x_analytical = [0, x_P, L];


phi_analytical = AnalyticalSolution('U', U, 'rho', rho, 'Gamma', Gamma, ...
                                    'L', L, 'x', x_analytical, ...
                                    'phi_A', phi_A, 'phi_B', phi_B);

phi_FVM = [phi_A;phi];
figure;
hold on;
plot(x_analytical, phi_FVM, 'rs-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', 'FVM (power-law differencing scheme)Numerical Solution');

plot(x_analytical, phi_analytical, 'b--', 'LineWidth', 2,'DisplayName', 'Analytical Solution');
title('1D Convection-Diffusion Solution (Power-Law Scheme)');
xlabel('Distance (m)', 'FontSize', 14);
ylabel('\phi', 'FontSize', 14);
title(sprintf('Numerical vs Analytical Solution (U = %.1f, N = %d)', U, N), 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;
% Visualize results
% plot(linspace(0,L,N+1), phi, '-o');
% title('1D Convection-Diffusion Solution (Power-Law Scheme)');
% xlabel('Position');
% ylabel('\phi');
% grid on;
