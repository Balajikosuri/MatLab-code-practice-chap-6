clc; clear; close all;

% Define parameters
h = 1;   % Distance between plates
U = 1;   % Velocity of upper plate
y = linspace(0, h, 10);  % Discretized y-axis values
u = (U / h) * y;  % Velocity profile equation u(y) = (U/h) * y

% Plot the stationary and moving plates
figure;
hold on;
plot([-0.2, 1.2], [0, 0], 'r', 'LineWidth', 3);  % Bottom plate (y = 0)
plot([-0.2, 1.2], [h, h], 'g', 'LineWidth', 3);  % Top plate (y = h)

% Plot vertical boundary lines (y = 0 and y = h)
plot([0, 0], [0, h], 'k', 'LineWidth', 2);  % Vertical line at u = 0
%plot([U, U], [0, h], 'k', 'LineWidth', 2);  % Vertical line at u = U

% Plot the velocity profile (linear distribution)
plot(u, y, 'b', 'LineWidth', 2); % Velocity profile u(y) vs y

% Add arrows starting from y=0 and y=h pointing to velocity profile
for i = 1:length(y)-1
    quiver(0, y(i), u(i), 0, 'LineWidth', 1.5, 'Color', 'k', 'MaxHeadSize', 0.2, 'AutoScale', 'off');

end

% Labels and title
xlabel('Velocity u(y)');
ylabel('Height y');
title('Plane Couette Flow - Velocity Profile');
axis([-0.2, U+0.2, -0.1, h+0.1]); % Set axis limits
grid on;
hold off;
