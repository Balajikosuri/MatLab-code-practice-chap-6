clear; clc; close all;

% Domain and grid
Lx = 3; Ly = 3;
Nx = 3; Ny = 3;
dx = Lx / Nx; dy = Ly / Ny;

% Boundary conditions
Phi_Left = 100; Phi_Right = 0;
Phi_Top = 0;    Phi_Bottom = 100;

% Physical properties
rho = 1; Gamma = 1;
u = 1; v = 4; % velocities

% Source terms (set to zero for simplicity)
S_u = 0; S_P = 0;

% QUICK requires 2 ghost cells on each side for a 3-point stencil
phi = zeros(Ny+4, Nx+4);

% Set boundary conditions (Dirichlet)
phi(3:Ny+2,1) = Phi_Left;       % West
phi(3:Ny+2,Nx+4) = Phi_Right;   % East
phi(1,3:Nx+2) = Phi_Top;        % North
phi(Ny+4,3:Nx+2) = Phi_Bottom;  % South

% Initialize ghost cells (mirror for Dirichlet)
phi(3:Ny+2,2) = Phi_Left;
phi(3:Ny+2,Nx+3) = Phi_Right;
phi(2,3:Nx+2) = Phi_Top;
phi(Ny+3,3:Nx+2) = Phi_Bottom;

% Iterative solver parameters
tol = 1e-6;
maxIter = 10000;
error = 1;
iter = 0;

% Precompute coefficients
Fe = rho * u * dy;
Fw = rho * u * dy;
Fn = rho * v * dx;
Fs = rho * v * dx;
De = Gamma * dy / dx;
Dw = Gamma * dy / dx;
Dn = Gamma * dx / dy;
Ds = Gamma * dx / dy;

while error > tol 
    phi_old = phi;
    for j = 3:Ny+2
        for i = 3:Nx+2
            % Skip boundary cells (handled by BCs)
            if i == 3 || i == Nx+2 || j == 3 || j == Ny+2
                % Use upwind for boundary-adjacent cells
                aE = De - max(0, Fe);
                aW = Dw + max(0, Fw);
                aN = Dn - max(0, Fn);
                aS = Ds + max(0, Fs);
                aP = aE + aW + aN + aS - S_P;
                phi(j,i) = (aE*phi(j,i+1) + aW*phi(j,i-1) + ...
                            aN*phi(j-1,i) + aS*phi(j+1,i) + S_u) / aP;
            else
                % Interior: QUICK scheme for convection terms
                % East face
                if u > 0
                    phi_e = (3/8)*phi(j,i+1) + (6/8)*phi(j,i) - (1/8)*phi(j,i-1);
                else
                    phi_e = (3/8)*phi(j,i) + (6/8)*phi(j,i+1) - (1/8)*phi(j,i+2);
                end
                % West face
                if u > 0
                    phi_w = (3/8)*phi(j,i-1) + (6/8)*phi(j,i-2) - (1/8)*phi(j,i-3);
                else
                    phi_w = (3/8)*phi(j,i-2) + (6/8)*phi(j,i-1) - (1/8)*phi(j,i);
                end
                % North face
                if v > 0
                    phi_n = (3/8)*phi(j-1,i) + (6/8)*phi(j-2,i) - (1/8)*phi(j-3,i);
                else
                    phi_n = (3/8)*phi(j-2,i) + (6/8)*phi(j-1,i) - (1/8)*phi(j,i);
                end
                % South face
                if v > 0
                    phi_s = (3/8)*phi(j+1,i) + (6/8)*phi(j,i) - (1/8)*phi(j-1,i);
                else
                    phi_s = (3/8)*phi(j,i) + (6/8)*phi(j+1,i) - (1/8)*phi(j+2,i);
                end

                % Discretized equation
                aE = De; aW = Dw; aN = Dn; aS = Ds;
                aP = aE + aW + aN + aS - S_P;

                phi(j,i) = (aE*phi_e + aW*phi_w + aN*phi_n + aS*phi_s + S_u) / aP;
            end
        end
    end
    error = max(max(abs(phi - phi_old)));
    iter = iter + 1;
end

disp(['Converged in ', num2str(iter), ' iterations. Final error: ', num2str(error)]);

% Extract internal field (remove ghost cells)
phi_internal = phi(3:Ny+2, 3:Nx+2);

% Plot results
x = linspace(dx/2, Lx-dx/2, Nx);
y = linspace(dy/2, Ly-dy/2, Ny);
[X, Y] = meshgrid(x, y);

figure;
contourf(X, Y, phi_internal, 20, 'LineColor', 'none');
colorbar; title('2D QUICK Scheme Solution');
xlabel('x'); ylabel('y');

figure;
surf(X, Y, phi_internal, 'EdgeColor', 'none');
colorbar; title('3D QUICK Scheme Solution');
xlabel('x'); ylabel('y'); zlabel('\phi');
view(45,30);

disp('Final phi (internal nodes):');
disp(phi_internal);
