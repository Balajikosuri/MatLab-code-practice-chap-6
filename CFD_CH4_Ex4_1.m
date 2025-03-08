clc; clear; close all;
% Temp distribution in 1D steady sate source free 
% Temperature Distribution: Numerical vs Analytical
% 1D steady-state heat conduction in a rod without internal heat generation

% Given Data
L = 0.5;       % Length of the rod (m)
N = 5;         % Number of internal nodes
T_A = 100;     % Left boundary temperature (°C)
T_B = 500;     % Right boundary temperature (°C)
k = 1000;      % Thermal conductivity (W/m·K)
A_cross = 0.01; % Cross-sectional area (m²)
S_p = 0;
% Compute grid spacing using dx = L/N
dx = L / N;

% Compute coefficients
a_W = (k * A_cross) / dx;  
a_E = (k * A_cross) / dx;  

% Initialize Coefficient Matrix A and RHS Vector B

[A,B] = TriDiagonlCoeffMatrix(N, a_W, a_E, T_A, T_B,S_p);

% Solve for Internal Node Temperatures (T)
T_internal = A \ B;

% Compute x-values (cell centers) using x_P = (P - 0.5) * dx
P = 1:N;  % Node indices
x_P = (P - 0.5) * dx;  % Apply correct positioning formula

% Include boundary points
x_full = [0, x_P, L];  
T_full = [T_A; T_internal; T_B];  % Full temperature distribution

% Compute Analytical Solution: T(x) = 800x + 100
% x_analytical = linspace(0, L, 100);
T_analytical = ((T_B-T_A)/L).*(x_full) + T_A;

disp('x_full')
disp(x_full)

% Display Matrices
fprintf('Coefficient Matrix A:\n');
disp(A);
fprintf('RHS Vector B:\n');
disp(B);
fprintf('Computed Internal Temperatures:\n');
disp(T_internal);

% Plot Numerical vs Analytical Solution
figure;
plot(x_full, T_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot(x_full, T_analytical, 'b-', 'LineWidth', 2);
hold off;

xlabel('Distance along the rod (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution: Numerical vs Analytical','FontSize', 14);
subtitle('1D steady-state heat conduction in a rod without internal heat generation', 'FontSize', 12);
grid on;
legend('Numerical Solution', 'Analytical Solution');

% Mark Internal Node Values on the Plot
for i = 1:N
    text(x_full(i+1), T_internal(i) + 10, sprintf('T%d = %.2f°C', i, T_internal(i)), 'FontSize', 10, 'Color', 'black');
end
