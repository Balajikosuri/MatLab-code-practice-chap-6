

L = 1; % Plate thickness (m)
T_inf = 20; % Ambient temperature (°C)
T_A = 100; % Boundary temperature at x=0
T_B = 0; % Boundary temperature at x=L
n_square = 25; % Thermal parameter (1/m²)
n = sqrt(n_square); % Derived parameter
x_analytical = 1;
T_analytical = T_inf + (T_A - T_inf) * (cosh(n * (L - x_analytical)) / cosh(n * L));
disp(T_analytical)