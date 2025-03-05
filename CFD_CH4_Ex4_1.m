% MATLAB Code to Solve Equation 4.23 (Steady-State 1D Heat Conduction)
clc; clear; close all;

% Define the coefficient matrix A
A = [-300  100   0    0    0;
      100 -200  100    0    0;
       0   100 -200  100    0;
       0    0   100 -200  100;
       0    0    0   100 -300];

% Define the right-hand side (RHS) vector b
b = [-200*100; 0; 0; 0; -200*500];

% Solve for Temperature Vector T using MATLAB's built-in solver
T_internal = A \ b;

% Display computed internal temperatures
fprintf('Computed Internal Temperatures:\n');
disp(T_internal);

% Define total number of nodes (including boundaries)
N = length(T_internal);  
L = 0.5;  % Length of the rod (m)
dx = L / N;  % Uniform spacing

% Compute x-values (cell centers) using the correct formula
P = 1:N;  % Node indices
x_P = 0 + (P - 0.5) * dx;  % Apply x_P = x_0 + (P - 1/2) * dx

% Include boundary points
x_full = [0, x_P, L];  
disp(['Number of x_full points: ', num2str(length(x_full))]);

% Full temperature distribution (including boundary temperatures)
T_full = [100; T_internal; 500];

% Analytical solution T(x) = 800x + 100
x_analytical = linspace(0, L, 100);
T_analytical = 800 * x_analytical + 100;
disp('x_analytical')
disp(x_analytical)

% Plot Temperature Distribution
figure;
plot(x_full, T_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot(x_analytical, T_analytical, 'b-', 'LineWidth', 2);
hold off;

xlabel('Distance along the rod (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution (Numerical vs Analytical)');
grid on;
legend('Numerical Solution', 'Analytical Solution');

% Mark Internal Node Values on the Plot
for i = 1:N
    text(x_full(i+1), T_internal(i) + 10, sprintf('T%d = %.2f°C', i, T_internal(i)), 'FontSize', 10, 'Color', 'black');
end
