clc; clear; close all;

% Define parameters
h = 1;                      % Normalized plate spacing (h = 1)
y = linspace(0, h, 100); 
% disp(y);

% Normalized height (y/h)

% Eckert-Prandtl number values
Ec_Pr_values = [0, 2, 4, 6, 8];

% Compute temperature profiles for different Ec*Pr values
figure;
hold on;
colors = lines(length(Ec_Pr_values)); % Generate distinct colors

for i = 1:length(Ec_Pr_values)
    Ec_Pr = Ec_Pr_values(i);
    
    % Corrected equation with element-wise multiplication (.*)
    T_norm = (y/h) + (Ec_Pr/2) .* (y/h) .* (1 - (y/h));  
    % where T_norm  = (T - T_0) / (T_1 - T_0)
  
    
    plot(T_norm, y/h, 'Color', colors(i, :), 'LineWidth', 2, 'DisplayName', sprintf('Ec \\cdot Pr = %d', Ec_Pr));
end

% Labels and formatting
xlabel('(T - T_0) / (T_1 - T_0)', 'FontSize', 12);
ylabel('y/h', 'FontSize', 12);
title('Temperature Distribution in Plane Couette Flow','Case I. When the plates are kept at different termperatures.', 'FontSize', 14);
legend('Location', 'best');
grid on;
set(gca, 'FontSize', 12);
hold off;
