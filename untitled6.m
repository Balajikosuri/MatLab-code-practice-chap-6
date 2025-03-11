% Given parameters
q = 1000;  % Heat generation (W/m^3)
k = 0.5;   % Thermal conductivity (W/m·K)

T_A = 100; % Temperature at x=0
T_B = 200; % Temperature at x=L
L = 0.02;  % Plate thickness (m)
N = 5;  

dx = L / N;  % Uniform spacing
disp(dx)
P = 1:N;
x_P = 0 + (P - 0.5) * dx;  % Apply x_P = x_0 + (P - 1/2) * dx
% disp(x_P)
% Include boundary points
x_analytical = [0, x_P, L];  


% Compute analytical temperature using the equation
T_analytical = (-q / (2 * k)) * x_analytical.^2 + ((T_B - T_A + (q / (2 * k)) * L^2) / L) * x_analytical + T_A;
plot(x_analytical,T_analytical)