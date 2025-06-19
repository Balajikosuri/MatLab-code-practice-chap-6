% Couette Flow using FTCS and Crank-Nicolson
clc; clear; close all

% Parameters
nu = 0.000217;          % Kinematic viscosity (m^2/s)
h = 0.04;               % Height (m)
N = 40;                 % Number of spatial nodes
dy = h / N;
y = linspace(0, h, N+1);
u0 = 40;                % Velocity at y = 0

% Time settings
dt_values = [0.002, 0.00232];
t_plot = [0.2, 0.4, 0.6, 1.0];

for k = 1:length(dt_values)
    dt = dt_values(k);
    d = nu * dt / dy^2;
    nt = ceil(max(t_plot)/dt) + 1;

    %% FTCS Method
    u_ftcs = zeros(N+1, nt);
    u_ftcs(1, :) = u0;
    for n = 1:nt-1
        for i = 2:N
            u_ftcs(i, n+1) = u_ftcs(i, n) + d * (u_ftcs(i+1, n) - 2*u_ftcs(i, n) + u_ftcs(i-1, n));
        end
        u_ftcs(1, n+1) = u0; % Bottom wall
        u_ftcs(end, n+1) = 0; % Top wall
    end

    %% Crank-Nicolson Method
    A = diag((1 + d) * ones(N-1, 1)) + ...
        diag((-d/2) * ones(N-2, 1), 1) + ...
        diag((-d/2) * ones(N-2, 1), -1);
    u_cn = zeros(N+1, nt);
    u_cn(1, :) = u0;
    
    for n = 1:nt-1
        b = zeros(N-1,1);
        for i = 2:N
            b(i-1) = (1 - d) * u_cn(i, n) + (d/2) * (u_cn(i+1, n) + u_cn(i-1, n));
        end
        b(1) = b(1) + (d/2)*u0;  % Bottom BC
        u_inner = A \ b;
        u_cn(2:N, n+1) = u_inner;
        u_cn(1, n+1) = u0;
        u_cn(end, n+1) = 0;
    end

    %% Plot FTCS
    figure
    for j = 1:length(t_plot)
        idx = round(t_plot(j)/dt) + 1;
        plot(u_ftcs(:, idx), y, 'DisplayName', ['t = ' num2str(t_plot(j))])
        hold on
    end
    title(['FTCS: \Delta t = ' num2str(dt) ', d = ' num2str(d)])
    xlabel('u (m/s)'); ylabel('y (m)')
    legend; grid on

    %% Plot Crank-Nicolson
    figure
    for j = 1:length(t_plot)
        idx = round(t_plot(j)/dt) + 1;
        plot(u_cn(:, idx), y, 'DisplayName', ['t = ' num2str(t_plot(j))])
        hold on
    end
    title(['Crank-Nicolson: \Delta t = ' num2str(dt) ', d = ' num2str(d)])
    xlabel('u (m/s)'); ylabel('y (m)')
    legend; grid on
end
