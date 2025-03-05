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
L = 0.5; % Length of the rod (m)

% Define x-coordinates including boundary points
x_full = linspace(0, L, N + 2); % Including boundaries at x=0 and x=L

% Full temperature distribution (including boundary temperatures)
T_full = [100; T_internal; 500];

% Plot Temperature Distribution
figure;
plot(x_full, T_full, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Distance along the rod (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution (Solution of Equation 4.23)');
grid on;
legend('Temperature Profile');

% Mark Internal Node Values on the Plot
for i = 1:N
    text(x_full(i+1), T_internal(i) + 10, sprintf('T%d = %.2f°C', i, T_internal(i)), 'FontSize', 10, 'Color', 'black');
end
