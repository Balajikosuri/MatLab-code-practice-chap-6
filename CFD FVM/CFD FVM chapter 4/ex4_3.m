 % Function 1: Analytical Solution
 function plot_analytical_solution(L, T_inf, T_A, n, color)
    % Function to plot the analytical temperature distribution
    % Inputs:
    % L - Plate thickness (m)
    % T_inf - Ambient temperature (°C)
    % T_A - Boundary temperature at x=0
    % n - Thermal parameter sqrt(n_square)
    % color - Line color for the plot (string)

    % Define x-coordinates for smooth analytical curve
    x_analytical = linspace(0, L, 100); 

    % Compute Analytical Temperature
    T_analytical = T_inf + (T_A - T_inf) * (cosh(n * (L - x_analytical)) / cosh(n * L));

    % Plot Analytical Solution
    plot(x_analytical, T_analytical, color, 'LineWidth', 2);
    hold on;
end
% Function 2: Numerical Solution (N = 5 and N = 10)
function plot_numerical_solution(N, L, T_inf, T_A, T_B, n_square, color, marker)
    % Function to compute and plot numerical solution using FVM
    % Inputs:
    % N - Number of nodes
    % L - Plate thickness (m)
    % T_inf - Ambient temperature (°C)
    % T_A - Boundary temperature at x=0
    % T_B - Boundary temperature at x=L
    % n_square - Thermal parameter
    % color - Line color
    % marker - Marker style for the plot
    
    % Compute Grid Spacing
    dx = L / N;
    P = 1:N;
    x_P = (P - 0.5) * dx;
    x_full = [0, x_P, L]; % Include boundary points

    % Compute FVM Coefficients
    a_W = 1/dx;  
    a_E = 1/dx;  
    S_p = n_square * T_inf * dx; 
    S_u = n_square * dx; 
    q_R = 0;

    % Solve Using Tridiagonal Matrix Algorithm
    [A, B] = TriDiagonlCoeffMatrix(N, a_W, a_E, T_A, T_B, S_p, S_u, q_R);
    T_internal = A \ B;
    
    % Full temperature distribution including boundaries
    T_B = T_internal(end);
    T_full = [T_A; T_internal; T_B];  

    % Plot Numerical Solution
    plot(x_full, T_full, marker, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', color);
    hold on;
end

clc; clear; close all;

% Define Parameters
L = 1; % Plate thickness (m)
T_inf = 20; % Ambient temperature (°C)
T_A = 100; % Boundary temperature at x=0
T_B = 0; % Boundary temperature at x=L
n_square = 25; % Thermal parameter (1/m²)
n = sqrt(n_square); % Derived parameter

% Plot Analytical Solution
plot_analytical_solution(L, T_inf, T_A, n, 'b-');

% Plot Numerical Solution for N = 5 (Coarse Grid)
plot_numerical_solution(5, L, T_inf, T_A, T_B, n_square, 'r', 'ro-');

% Plot Numerical Solution for N = 10 (Fine Grid)
plot_numerical_solution(10, L, T_inf, T_A, T_B, n_square, 'g', 'gs-');

% Final Plot Formatting
title('Temperature Distribution: Analytical vs Numerical (FVM)', 'FontSize', 14);
xlabel('Distance (m)');
ylabel('Temperature (°C)');
legend('Analytical Solution', 'Numerical (FVM, N=5)', 'Numerical (FVM, N=10)');
grid on;
hold off;
