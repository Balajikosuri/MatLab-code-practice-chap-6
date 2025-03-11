clc; clear; close all;
%  1 D steady state Diffusion equation problem with sources
% Temperature Distribution in a Rod with Internal Heat Generation (Example 4.2)


% Define x-coordinates (cell centers)
L = 1; % Plate thickness (m)
N = 5; % Number of nodes
n_square = 25; %(/m^2 kA is constant and 
% Compute grid spacing using dx = L/N
dx = L / N;

% temperatures
T_inf = 20; %Ambient temperature 
T_A = 100; % Boundary temperature at x=0
T_B = 0; % Boundary temperature at x=L

a_W = 1/dx;  
a_E = 1/dx;  
S_p = n_square*T_inf*dx; 
S_u = n_square*dx; 
q_R = 0;

[A, B] = TriDiagonlCoeffMatrix(N, a_W, a_E, T_A, T_B,S_p,S_u,q_R);
disp("Matrix A")
disp(A)
disp("Matrix B")
disp(B)

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

c_1 = (T_B - T_A)/L + (q*L)/(2*k);
c_2 = T_A;
% disp('c2')
% disp(c_2)
% Compute analytical temperature usin*g the equation
T_analytical = -q*(x_analytical.^2)/(2*k) + c_1*x_analytical+c_2;
disp("T_analytical")
disp(T_analytical)


% Plot results
figure;
plot(x_analytical*100, T_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold on;
plot(x_analytical*100, T_analytical, 'b-', 'LineWidth', 2);
hold off;
title('Temperature Distribution: Numerical vs Analytical','FontSize', 14);
subtitle('Temperature Distribution in a Rod with Internal Heat Generation (Example 4.2)', 'FontSize', 12);

xlabel('Distance (cm)');
ylabel('Temperature (°C)');

legend('Numerical (FVM)', 'Analytical');
grid on;

% Annotate values
% for i = 1:N
%     text(x_analytical(i+1), T_internal(i) + 5, sprintf('T%d = %.2f°C', i, T_internal(i)), 'FontSize', 10, 'Color', 'black');
% end
