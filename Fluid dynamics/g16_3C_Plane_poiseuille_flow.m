clc; clear; close all;

% Define parameters
h = 2;   % Distance between plates
U = 0;   % Velocity of upper plate
mu = 4;
P = -4;
y = linspace(0, h/2, 10);  % Discretized y-axis values
u = -((P*h^2)/(8*mu))*(1-4*(y/h).^2);  % Velocity profile equation u(y) = (U/h) * y

% Plot the stationary and moving plates
figure;
hold on;
plot([-0.2, 1], [-h/2, -h/2], 'r', 'LineWidth', 3);  % Bottom plate (y = -h/2)
plot([-0.2, 1], [h/2, h/2], 'r', 'LineWidth', 3);  % Top plate (y = h/2)

% Plot vertical boundary lines (y = 0 and y = h)
plot([0, 0], [0, h/2], 'k', 'LineWidth', 2);  % Vertical line at u = 0
plot([0, 0], [0, -h/2], 'k', 'LineWidth', 2);  % Vertical line at u = 0

plot([0, 1], [0, 0], 'k', 'LineWidth', 2);  
% Add arrow at the end of the x-axis using quiver
quiver(1, 0, 0.05, 0, 'MaxHeadSize', 0.5, 'LineWidth', 2, 'Color', 'k'); 
% Label for X-axis
text(1, 0, '  X-Axis', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');



% Plot the velocity profile (linear distribution)
plot(u, y, 'b', 'LineWidth', 2); % Velocity profile u(y) vs y
plot(u, -y, 'b', 'LineWidth', 2); % Velocity profile u(y) vs y


% Add arrows starting from y=0 and y=h pointing to velocity profile
% for i = 1:length(y)-1
%     quiver(0, y(i), u(i), 0, 'LineWidth', 1.5, 'Color', 'k', 'MaxHeadSize', 0.2, 'AutoScale', 'off');
% 
% end

% Labels and title
xlabel('Velocity u(y)');
ylabel('Height y');
title('Plane poiseuille flow - Velocity Profile');
axis([-0.2, 2, -0.2, 2]); % Set axis limits
grid on;
hold off;
