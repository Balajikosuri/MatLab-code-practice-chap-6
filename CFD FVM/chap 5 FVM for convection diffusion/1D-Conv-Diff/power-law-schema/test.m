clc; clear;

% Parameters
L = 1;           % Domain length
N = 5;           % Number of control volumes
rho = 1;         % Density
u = 1;           % Velocity
Gamma = 0.1;     % Diffusion coefficient
phi_0 = 1;       % Left boundary value
phi_L = 0;       % Right boundary value

dx = L / N;                        % Uniform cell width
x = linspace(dx/2, L-dx/2, N);     % Cell centers
F = rho * u;                       % Convection flux (assumed constant)

% Initialize arrays
aP = zeros(1, N);
aE = zeros(1, N);
aW = zeros(1, N);
Su = zeros(1, N);

% Loop over internal nodes (1 to N)
for i = 1:N
    % Diffusion conductance varies at boundaries
    if i == 1
        De = Gamma / dx;          % East face (node 1 to 2)
        Dw = Gamma / (dx / 2);    % West face (boundary to node 1)
    elseif i == N
        De = Gamma / (dx / 2);    % East face (node N to boundary)
        Dw = Gamma / dx;          % West face (node N-1 to N)
    else
        De = Gamma / dx;          % East face
        Dw = Gamma / dx;          % West face
    end

    % Péclet numbers at faces
    Pe_e = F / De;
    Pe_w = F / Dw;

    % Power-law factors
    fE = max([0, (1 - 0.1 * abs(Pe_e))^5]);
    fW = max([0, (1 - 0.1 * abs(Pe_w))^5]);

    aE(i) = De * fE;
    aW(i) = Dw * fW;

    % Central node coefficient
    aP(i) = aE(i) + aW(i);

    % Boundary contributions
    if i == 1
        Su(i) = aW(i) * phi_0;
    elseif i == N
        Su(i) = aE(i) * phi_L;
    end
end

% Construct system matrix A and RHS vector b
A = zeros(N, N);
b = Su';

for i = 1:N
    if i ~= 1
        A(i, i-1) = aW(i);
    end
    A(i, i) = aP(i);
    if i ~= N
        A(i, i+1) = aE(i);
    end
end

% Solve the system
phi = A \ b;

% Add boundary values for plotting
x_full = [0, x, L];
phi_full = [phi_0, phi', phi_L];

% Plotting
plot(x_full, phi_full, '-o', 'LineWidth', 2);
xlabel('x'); ylabel('\phi');
title('1D Steady Convection-Diffusion using Power-Law Scheme (Patankar)');
grid on;
