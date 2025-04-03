% Simulating Fluid Flow in MATLAB with 3D Visualization

% Number of grid points
N = 3+2;

% Domain size
L = 1;

% Corrected grid spacing
h = L / N; 

% Iteration number
iterations = 0;

% Thermal conductivity
k = 1;

% Cross-section area
A = 1;

% Initializing temperature field
T = zeros(N, N);


% Setting boundary conditions
T(:, 1) = 400; % Left boundary TL
T(:, end) = 300; % Right boundary TR
T(1, :) = 500; % Top boundary TT
T(end, :) = 200; % Bottom boundary TB
disp("Initial T ")
disp(T)
% Initializing iterated temperature
T_new = T;

% Error related variables
epsilon = 1.E-8;
numerical_error = 1;

% Plot for numerical error
figure(10);

% Checking the error tolerance
while numerical_error > epsilon
    % Computing for all interior points
    for i = 2:N-1
        for j = 2:N-1
            % fprintf('The value of (i,j)= (%d,%d)\n',i,j)
            a_E = k * A / h;
            a_W = k * A / h;
            a_N = k * A / h;
            a_S = k * A / h;
            a_P = a_E + a_W + a_N + a_S;
            T_new(i, j) = (a_E * T(i, j+1) + a_W * T(i, j-1) + a_N * T(i-1, j) + a_S * T(i+1, j)) / a_P;
        end
    end
    
    % Resetting the numerical error and recalculate
    numerical_error = sum(sum(abs(T(2:N-1, 2:N-1) - T_new(2:N-1, 2:N-1))));
    
    % Iteration advancement and reassignment
    iterations = iterations + 1;
    T = T_new;
    
    % Plotting numerical error
    if mod(iterations, 1000) == 0
        figure(10);
        semilogy(iterations, numerical_error, 'ko');
        pause(0.01);
    end
end
disp("T_new")
disp(T_new)
% 3D Surface Plot of Temperature Distribution
disp('---------------------')
% Defining the position vector and the grid from indices
x_dom = linspace(0, L, N);
y_dom = linspace(0, L, N);
[X, Y] = meshgrid(x_dom, y_dom);

% Plotting the 3D surface
figure(11);
contourf(X, Y, T, 12);

surf(X, Y, T, 'EdgeColor', 'none'); % Smooth 3D surface plot
colormap jet; % Use jet colormap for better visualization
colorbar; % Add color scale
xlabel('X-Domain');
ylabel('Y-Domain');
zlabel('Temperature T(x,y)');
title('3D Temperature Distribution');
view(45, 30); % Adjust viewing angle for better visibility
grid on;