clc; clear; close all;

% Define parameters
h = 1;                      % Normalized plate spacing (h = 1)
y = linspace(-h, h, 100);   % Discretized y-values from -h to h (centered at y = 0)
T_norm = 1 - (16 .* (y.^4 / h.^4));  % Dimensionless temperature profile

% Corrected fprintf statement
% fprintf('T_norm values (first 5): %.4f, %.4f, %.4f, %.4f, %.4f\n', T_norm(1:5));
% fprintf('y values (first 5): %.4f, %.4f, %.4f, %.4f, %.4f\n', y(1:5));

% Plot the temperature distribution
figure;

hold on;
plot(T_norm, y/(2*h), 'k', 'LineWidth', 2); % Corrected normalization

% X-axis line
plot([0, 1], [0, 0], 'k', 'LineWidth', 2); 

% Markers for each 0.2 unit on x-axis
x_marks = 0:0.2:1;  % Mark positions at 0, 0.2, 0.4, ..., 1
y_marks = zeros(size(x_marks));  % All y-values at zero for markers
plot(x_marks, y_marks, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % Black circular markers
for i = 1:length(x_marks)
    text(x_marks(i), -0.05, sprintf('%.1f', x_marks(i)), 'FontSize', 12, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end
text(0.5, 0, 'T_m', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% 
xlabel('$(T - T_o)/(T_m - T_o)$', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('$y/2h$', 'FontSize', 14, 'Interpreter', 'latex');

% Title and formatting
title({'Dimensionless Temperature Distribution in Poiseuille Flow', ...
       'Steady Laminar Flow of an Incompressible Fluid'}, ...
       'FontSize', 14);

grid on;
axis([-0.2, 1.2, -0.6, 0.6]); % Adjusted axis for clarity

% Set x-axis ticks from 0 to 1 with 0.2 spacing
xticks(0:0.2:1);  

set(gca, 'FontSize', 12);

% Annotations for axes

annotation('textarrow', [0.1 0.9], [0.05 0.05], 'FontSize', 14);

annotation('textarrow', [0.1 0.9], [0.05 0.05], 'FontSize', 14);

hold off;
