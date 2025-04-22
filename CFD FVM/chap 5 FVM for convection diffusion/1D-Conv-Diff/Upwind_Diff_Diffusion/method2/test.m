% Parameters
L = 1.0;        % Length of the domain
N = 5;         % Number of control volumes
rho = 1.0;      % Density
u = 1.0;        % Velocity
Gamma = 0.1;    % Diffusion coefficient

phi_A = 1;      % Boundary condition at x=0
phi_B = 0;      % Boundary condition at x=L

dx = L / N;     % Grid spacing
F = rho * u;    % Convective flux
D = Gamma / dx; % Diffusive conductance

% Initialize matrices
A = zeros(N+1, N+1);
b = zeros(N+1, 1);

% Boundary nodes
A(1,1) = 1;
b(1) = phi_A;
A(end,end) = 1;
b(end) = phi_B;

% Internal nodes: Upwind scheme
for i = 2:N
    aW = D + max(F, 0);
    aE = D + max(-F, 0);
    aP = aW + aE;
    A(i, i-1) = -aW;
    A(i, i)   = aP;
    A(i, i+1) = -aE;
end

% Solve system
phi = A \ b;

%% 
P = 1:N;
x_ps = (P-1/2)*dx;
x = [0,x_ps,L];
plot(x, phi, '-o');
xlabel('x');
ylabel('\phi');
title('Steady-State 1D Convection-Diffusion (Upwind FVM)');
grid on;
