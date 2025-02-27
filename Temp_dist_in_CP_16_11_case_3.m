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

T_norm = 1-(y.^2/h.^2); 

% where T_norm = (T-T1)/(Tr-T1)
% y = h/2 
    
plot(T_norm, y/h, 'k', 'LineWidth', 2);
annotation('textarrow', [0.1, 0.1], [0.3, 0.6],'FontSize', 12, 'FontWeight', 'bold', 'HeadStyle', 'vback2','LineWidth', 2, 'Color', 'k');
annotation('textarrow', [0.3, 0.6], [0.1, 0.1], 'FontSize', 12, 'FontWeight', 'bold', 'HeadStyle', 'vback2', 'LineWidth', 2, 'Color', 'k');

% 


text(0, h/2, num2str(h/2), 'VerticalAlignment', 'bottom', 'FontSize', 12);
text(h/2,0, num2str(h/2), 'VerticalAlignment', 'bottom', 'FontSize', 12);


% Labels and formatting
xlabel('(T-T1)/(Tr-T1)', 'FontSize', 12);
ylabel('y/h', 'FontSize', 12);

title('Temperature Distribution in Plane Couette Flow',{'Case III: When at one of the plates,say the lower stationary plate, no heat transfer takes place (adiabatic wall),'
       'and the other moving wall is kept at temperature T1'}, ...
       'FontSize', 14);

% legend('Location', 'best');
grid on;
set(gca,'FontSize', 12);
axis([-0.2, 2, -0.2, 2]);
hold off;
