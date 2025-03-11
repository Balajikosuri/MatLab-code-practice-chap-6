% MATLAB Code to Solve Equation 4.23 (Steady-State 1D Heat Conduction)
clc; clear; close all;

% Given Data
L = 0.02; % Plate thickness (m)
k = 0.5; % Thermal conductivity (W/m.K)
q = 1000; % Heat generation (W/m^3)
N = 5; % Number of control volumes
dx = L / N; % Grid spacing

% Compute x-values (cell centers) using correct FVM centering formula
P = 1:N;  % Node indices
x_P = (P - 0.5) * dx; % Apply x_P = x_0 + (P - 1/2) * dx

% Include boundary points
x_full = [0, x_P, L];  
disp('Computed x_P values:');
disp(x_P);

% Define Coefficient Matrix A and RHS vector B
A = zeros(N, N);
B = zeros(N, 1);

% Interior node equations
for i = 2:N-1
    A(i, i-1) = 1; % West node
    A(i, i) = -2;  % Central node
    A(i, i+1) = 1; % East node
    B(i) = - (q * dx^2) / k; % Source term
end

% Boundary conditions
T_A = 100; % At x = 0
T_B = 200; % At x = L

% First node equation (including T_A)
A(1,1) = -2; A(1,2) = 1;
B(1) = - (q * dx^2) / k - T_A; 

% Last node equation (including T_B)
A(N,N) = -2; A(N,N-1) = 1;
B(N) = - (q * dx^2) / k - T_B;

% Solve for internal temperatures
T_internal = A \ B;
disp(A)

% Full temperature distribution including boundaries
T_full = [T_A; T_internal; T_B];

% Analytical Solution
x_analytical = linspace(0, L, 100);
T_analytical = -1e6 * x_analytical.^2 + 25000 * x_analytical + 100;

% Plot results
figure;
plot(x_full, T_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot(x_analytical, T_analytical, 'b-', 'LineWidth', 2);
hold off;

xlabel('Position (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution: Numerical vs Analytical');
legend('Numerical (FVM)', 'Analytical');
grid on;

% Annotate values
for i = 1:N
    text(x_full(i+1), T_internal(i) + 5, sprintf('T%d = %.2f°C', i, T_internal(i)), 'FontSize', 10, 'Color', 'black');
end
