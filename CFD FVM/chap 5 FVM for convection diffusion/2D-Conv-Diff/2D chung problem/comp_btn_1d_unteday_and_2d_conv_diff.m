clc; clear; close all

%% Parameters for both problems
nu = 0.000217;      % Kinematic viscosity
h = 0.04;           % Height of channel
u0 = 40;            % Bottom wall velocity (Couette) / Inlet BC (convection)
Ny = 41;            % Number of grid points in y
dy = h / (Ny - 1);
y = linspace(0, h, Ny);

t_plot = [0.2, 0.4, 0.6, 1.0];   % time/x locations to sample

%% -------------------------------
%% Part 1: 1D Couette Flow (FTCS)
%% -------------------------------
dt = 0.002;
Nt = ceil(max(t_plot) / dt) + 1;
d = nu * dt / dy^2;

u_ftcs = zeros(Ny, Nt);
u_ftcs(1, :) = u0;   % Bottom wall BC

for n = 1:Nt-1
    for j = 2:Ny-1
        u_ftcs(j, n+1) = u_ftcs(j, n) + d * (u_ftcs(j+1, n) - 2*u_ftcs(j, n) + u_ftcs(j-1, n));
    end
    u_ftcs(1, n+1) = u0;
    u_ftcs(end, n+1) = 0;
end

%% ---------------------------------------
%% Part 2: 2D Convection-Diffusion (Steady)
%% ---------------------------------------
L = 1.0;                 % Length in x-direction
Nx = 1000;
dx = L / (Nx - 1);
x = linspace(0, L, Nx);
U = 1.0;                 % Constant x-velocity

d = nu * dx / (U * dy^2);
fprintf('The Diffusion Number/ The Numerical Stability (d) = %.4f\n', d)
if d <1/2
    fprintf('The solution is Stable (i.e d = %.4f )for Nx = %.4f\n', d,Nx);
else
    fprintf('Unstabel solution(i.e d = %.4f )for Nx = %f ) = %.4f\n', d,Nx);
end

phi = zeros(Ny, Nx);
phi(1,1) = u0;           % φ = 40 at bottom inlet (analogous to initial u=40)
phi(2:end,1) = 0;        % Elsewhere inlet = 0

for i = 1:Nx-1
    for j = 2:Ny-1
        phi_yy = (phi(j+1,i) - 2*phi(j,i) + phi(j-1,i)) / dy^2;
        phi(j,i+1) = phi(j,i) + dx/U * nu * phi_yy;
    end
    phi(1,i+1) = u0;    % bottom wall
    phi(end,i+1) = 0;   % top wall
end

%% -------------------------
%% Comparison Plotting
%% -------------------------
x_sample = t_plot;  % Since x = t in analogy

figure
for k = 1:length(t_plot)
    t_idx = round(t_plot(k)/dt) + 1;
    x_idx = round(x_sample(k)/dx) + 1;

    subplot(2,2,k)
    plot(u_ftcs(:, t_idx), y, 'b-o', 'DisplayName', ['1D FTCS at t = ' num2str(t_plot(k))])
    hold on
    plot(phi(:, x_idx), y, 'r--s', 'DisplayName', ['2D Convection at x = ' num2str(x_sample(k))])
    xlabel('Velocity / \phi (m/s)')
    ylabel('y (m)')
    title(['Comparison at t = x = ' num2str(t_plot(k))])
    legend; grid on
end

sgtitle('Comparison: 1D Couette (FTCS) vs 2D Convection-Diffusion')
