clc; clear; close all;

% Define parameters
h = 1;                      % Normalized plate spacing (h = 1)
y = linspace(0, h, 100);
U = 2;
mu = 2;
k = 1;  % Thermal conductivity (W/m·K)

% h = 1;                      % Normalized plate spacing (h = 1)
% y = linspace(0, h, 100);    % Normalized height (y/h)
% U = 2;                      % Velocity of moving plate (m/s)
% mu = 0.001;                 % Dynamic viscosity (Pa·s)
% k = 0.6;                    % Thermal conductivity (W/m·K)
% T0 = 300;                   % Initial temperature (K)
% Tm = 350;


% Compute temperature profiles for different  values
figure;
hold on;
% (T-To) = (mu*U.^2/2*k).*(y/h).*(1-(y/h));

plot([0, h], [0, 0], 'k', 'LineWidth', 3); % Lower plate
plot([0, h], [h, h], 'k', 'LineWidth', 3); % Upper plate

text(1, h, 'Upper Plate (y = h, U = 0)', 'FontSize', 12, 'FontWeight', 'bold');
text(1, 0, 'Lower Plate (y = 0, U = 0)', 'FontSize', 12, 'FontWeight', 'bold');


% 
T_norm = 4*(y/h).*(1-(y/h)); 
% where T_norm = (T-To)/(Tm-To)
% y = h/2 

    
plot(T_norm, y/h, 'k', 'LineWidth', 2);
annotation('textarrow', [0.1, 0.1], [0.3, 0.6],'FontSize', 12, 'FontWeight', 'bold', 'HeadStyle', 'vback2','LineWidth', 2, 'Color', 'k');
annotation('textarrow', [0.3, 0.6], [0.1, 0.1], 'FontSize', 12, 'FontWeight', 'bold', 'HeadStyle', 'vback2', 'LineWidth', 2, 'Color', 'k');

% 

plot([0,1], [h/2,h/2], 'k', 'LineWidth', 2);
text(0, h/2, num2str(h/2), 'VerticalAlignment', 'bottom', 'FontSize', 12);
text(1, h/2, num2str(h/2), 'VerticalAlignment', 'bottom', 'FontSize', 12);
text(h/2, h/2, 'T_m', 'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

% Labels and formatting
xlabel('(T-To)/(Tm-To)', 'FontSize', 12);
ylabel('y/h', 'FontSize', 12);
title('Temperature Distribution in Plane Couette Flow','Case II.When both the plates are kept at the same constant temperature To', 'FontSize', 14);
% legend('Location', 'best');
grid on;
set(gca,'FontSize', 12);
axis([0, 2, -0.2, 2]);
hold off;
