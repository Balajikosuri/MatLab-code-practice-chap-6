clear; clc; close all;

% Domain size
Lx = 1.0;  % Length in x-direction (m)
Ly = 1.0;  % Length in y-direction (m)

% Grid parameters
Nx = 4;  % Number of control volumes in x-direction
Ny = 4;  % Number of control volumes in y-direction
dx = Lx / Nx-1;
dy = Ly / Ny-1;

% Thermal conductivity
k = 1000; % W/(m·K)

% Initialize temperature field
T = zeros(Ny+1, Nx+1);

% Apply Dirichlet Boundary Conditions
T(:, 1)   = 400; % Left boundary (x = 0)
T(:, end) = 300; % Right boundary (x = Lx)
T(1, :)   = 500; % Top boundary (y = Ly)
T(end, :) = 200; % Bottom boundary (y = 0)
fprintf('The Initial Matrix With all Boundaries(Top Left Right Bottom Left \n')
disp(T)
% Finite Volume Coefficients
aW = k;
aE = k;
aN = k ;
aS = k ;

% Modify coefficients near boundaries
aW_b = 2 * k; 
aE_b = 2 * k; 
aN_b = 2 * k ; 
aS_b = 2 * k ; 

% Iterative solver (Gauss-Seidel)
tol = 1e-6;
error = 1;
maxIter = 10000;
iter = 0;

while error > tol && iter < maxIter
    T_old = T;

    % Update interior points
    for i = 2:Nx
        for j = 2:Ny
            
            if i == 2  % Near left boundary
                aW_eff = aW_b;
                % fprintf('(%d,%d)',j,i)
                % fprintf('\n the aP = %d\n)',aP)
            else
                aW_eff = aW;
            end

            if i == Nx  % Near right boundary
                aE_eff = aE_b;
                % fprintf('(%d,%d)',j,i)
            else
                aE_eff = aE;
            end

            if j == 2  % Near bottom boundary (edit to top)
                % fprintf('(%d,%d)',j,i)
                aN_eff = aN_b;
            else
                aN_eff = aN;
            end

            if j == Ny  % Near top boundary (edit to bottom)
                fprintf('(%d,%d)',j,i)
                aS_eff = aS_b;
            else
                aS_eff = aS;
            end

            % Compute new temperature using FVM equation
            aP = aW_eff + aE_eff + aN_eff + aS_eff;
            T(j, i) = (aE_eff * T(j, i+1) + aW_eff * T(j, i-1) + ...
                       aN_eff * T(j+1, i) + aS_eff * T(j-1, i)) / aP;
        end
    end
    
    % Compute error
    error = max(max(abs(T - T_old)));
    iter = iter + 1;
end
fprintf('The Final T Matrix Temp at respective nodes internal + nearer boundary\n')
disp(T)
disp(['Converged in ', num2str(iter), ' iterations with error ', num2str(error)]);


% Till here function is fine function call has to writen the T as Phi ,
% please use phi insted of T

% Plot solution
figure;
contourf(flipud(T), 12, 'LineColor', 'none'); 
colorbar;
clim([200 500]); 
colormap(jet);
title('2D Heat Conduction - Finite Volume Method');
xlabel('x (m)');
ylabel('y (m)');
