function RunCentrealDiffSchConvectionDiffusion(varargin)
    % Parse Name-Value pairs
    p = inputParser;
    
    addParameter(p, 'U', 1);               % Velocity
    addParameter(p, 'N', 10);              % Number of control volumes
    addParameter(p, 'L', 1);               % Length of the domain
    addParameter(p, 'Gamma', 0.1);         % Diffusion coefficient
    addParameter(p, 'rho', 1);             % Density
    addParameter(p, 'phi_A', 1);           % Left boundary value
    addParameter(p, 'phi_B', 0);           % Right boundary value
    % addParameter(p, 'AnalyticalSolution', NaN); %this a function of x/y/z pass the respective x/y/z value

    parse(p, varargin{:});
    U = p.Results.U;
    N = p.Results.N;
    L = p.Results.L;
    Gamma = p.Results.Gamma;
    rho = p.Results.rho;
    phi_A = p.Results.phi_A;
    phi_B = p.Results.phi_B;
    % AnalyticalSolution = p.Results.AnalyticalSolution; % its a fucunction pass x/y/z value in it 

    % Discretization
    dx = L / N;
    P = 1:N;
    x_P = (P - 0.5) * dx;
    x_analytical = [0, x_P, L];

    % Generate Tridiagonal Matrix
    [A, B] = TriDiagonalCoeffMatrix('N', N, ...
                                    'Diffusion', Gamma / dx, ...
                                    'Convection', rho * U, ...
                                    'PhiLeft', phi_A, ...
                                    'PhiRight', phi_B);

    phi_numerical = A \ B;
    phi_FVM_Num_full = [phi_A; phi_numerical; phi_B];

    % Analytical solution
    phi_analytical = AnalyticalSolution('U', U, 'rho', rho, 'Gamma', Gamma, ...
                                    'L', L, 'x', x_analytical, ...
                                    'phi_A', phi_A, 'phi_B', phi_B);
    % Plot
    figure;
    hold on;
    plot(x_analytical, phi_FVM_Num_full, 'rs-', 'LineWidth', 2, ...
         'MarkerSize', 8, 'DisplayName', 'FVM Numerical Solution');
    plot(x_analytical, phi_analytical, 'b--', 'LineWidth', 2, ...
         'DisplayName', 'Analytical Solution');
    xlabel('Distance (m)', 'FontSize', 14);
    ylabel('\phi', 'FontSize', 14);
    title(sprintf('Numerical vs Analytical Solution (U = %.1f, N = %d)', U, N), 'FontSize', 14);
    legend('Location', 'best');
    grid on;
    hold off;

    % Error Table
    Difference = abs(phi_FVM_Num_full - phi_analytical');
    Percentage_Error = (Difference ./ phi_analytical') * 100;

    T = table((0:N+1)', x_analytical', phi_FVM_Num_full, phi_analytical', ...
              Difference, Percentage_Error, ...
              'VariableNames', {'Node', 'Distance_m', 'FVM_Solution', ...
                                'Analytical_Solution', 'Difference', 'Percentage_Error'});

    fprintf("\nThe No.Of Nodes: %d, Velocity: %.4f m/s\n", N, U)
    fprintf('\n%-6s %-12s %-20s %-20s %-15s %-15s\n', 'Node', 'Distance(m)', ...
            'FVM Solution', 'Analytical Solution', 'Difference', 'Percentage Error');
    fprintf(repmat('-', 1, 100)); fprintf('\n');

    for k = 1:height(T)
        fprintf('%-4d %-12.4f %-20.6f %-20.4f %-15.3f %-15.2f\n', ...
                T.Node(k), T.Distance_m(k), T.FVM_Solution(k), ...
                T.Analytical_Solution(k), T.Difference(k), T.Percentage_Error(k));
    end
    disp('------ End of the Table ------------')
end

function phi = AnalyticalSolution(varargin)
    % Parse Name-Value pairs
    p = inputParser;

    addParameter(p, 'U', 1);               % Velocity
    addParameter(p, 'rho', 1);            % Density
    addParameter(p, 'Gamma', 0.1);        % Diffusion coefficient
    addParameter(p, 'L', 1);              % Domain length
    addParameter(p, 'x', []);             % x locations where phi is evaluated
    addParameter(p, 'phi_A', 1);          % Left boundary value
    addParameter(p, 'phi_B', 0);          % Right boundary value

    parse(p, varargin{:});
    U      = p.Results.U;
    rho    = p.Results.rho;
    Gamma  = p.Results.Gamma;
    L      = p.Results.L;
    x      = p.Results.x;
    phi_A  = p.Results.phi_A;
    phi_B  = p.Results.phi_B;

    % Analytical formula for 1D steady convection-diffusion problem
    exponent = (rho * U * x) / Gamma;
    denominator = exp((rho * U * L) / Gamma) - 1;
    numerator = exp(exponent) - 1;

    phi = phi_A + (phi_B - phi_A) * (numerator / denominator);
end


% phi_analytical = AnalyticalSolution('U', U, 'rho', rho, 'Gamma', Gamma, ...
                                    % 'L', L, 'x', x_analytical, ...
                                    % 'phi_A', phi_A, 'phi_B', phi_B);


clc; clear; close all;

% Physical and boundary values
L = 1.0;
Gamma = 0.1;
rho = 1.0;
phi_A = 1;
phi_B = 0;

U = 0.1;
N = 5;
U_values = [0.1,2.5];
N_values = [5,20];

for i = 1:length(U_values)
    for j = 1:length(N_values)
        U = U_values(i); % Current velocity
        N = N_values(j); % Current number of nodes
        RunCentrealDiffSchConvectionDiffusion('U', U, 'N', N, 'L',L, 'Gamma', Gamma, 'rho', rho, 'phi_A', phi_A, 'phi_B', phi_B);
    end
end
