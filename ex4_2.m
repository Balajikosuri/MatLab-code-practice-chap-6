clc; clear; close all;
%  1 D steady state Diffusion equation problem with sources 

% Define x-coordinates (cell centers)
L = 0.02; % Plate thickness (m)
N = 5; % Number of nodes
k = 0.5; % Thermal conductivity (W/m.K)
q = 1000*1000; % Heat generation (W/m^3)
A = 1;
% Compute grid spacing using dx = L/N
dx = L / N;

T_A = 100; % Boundary temperature at x=0
T_B = 200; % Boundary temperature at x=L

a_W = (k*A)/dx;  
a_E = (k*A)/dx;  
S_P = q*A*dx; 

[A, B] = TriDiagonlCoeffMatrix(N, a_W, a_E, T_A, T_B,S_P);
% Solve for internal temperatures
T_internal = A \ B;

% Display computed internal temperatures
fprintf('Computed Internal Temperatures:\n');
disp(T_internal);



dx = L / N; % Grid spacing
P = 1:N; % Node indices
x_P = (P - 0.5) * dx;  % Apply x_P = x_0 + (P - 1/2) * dx

% Include boundary points
x_analytical = [0, x_P, L]; 
disp("x_analytical")
disp(x_analytical)

% Full temperature distribution including boundaries
T_full = [T_A; T_internal; T_B ];  

% Compute analytical temperature using the equation
T_analytical = -10.^6 * (x_analytical).^2+25000*x_analytical+100;
disp("T_analytical")
disp(T_analytical)

% Plot results
figure;
plot(x_analytical, T_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
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
    text(x_analytical(i+1), T_internal(i) + 5, sprintf('T%d = %.2f°C', i, T_internal(i)), 'FontSize', 10, 'Color', 'black');
end
