clc; clear; close all;

% Given Parameters
L = 1.0;        % Length of the domain
Gamma = 0.1;    % Diffusion coefficient (Γ)
rho = 1.0;      % Density

phi_A = 1; % Left boundary condition
phi_B = 0; % Right boundary condition

% Define different values of U (velocity) and N (grid points)
U_values = [0.1, 2.5];  % Cases: Low and high velocity
N_values = [5, 20];     % Cases: Coarse and fine grid

%function declaration : 

function SolveCDsConvectionDiffusionCase(varargin)
    % Parse input arguments using key-value pairs
    p = inputParser;
    addParameter(p, 'U', 0.1);
    addParameter(p, 'N', 5);
    addParameter(p, 'L', 1.0);
    addParameter(p, 'Gamma', 0.1);
    addParameter(p, 'rho', 1.0);
    addParameter(p, 'phi_A', 1);
    addParameter(p, 'phi_B', 0);
    addParameter(p, 'anaylitical_fun', NaN);

    parse(p, varargin{:});
    U = p.Results.U;
    N = p.Results.N;
    L = p.Results.L;
    Gamma = p.Results.Gamma;
    rho = p.Results.rho;
    phi_A = p.Results.phi_A;
    phi_B = p.Results.phi_B;
    anaylitical_fun = p.Results.anaylitical_fun;

    dx = L / N;

    % Cell-centered grid
    P = 1:N;
    x_P = (P - 0.5) * dx;
    x_analytical = [0, x_P, L];

    % Build matrix
    [A, B] = TriDiagonalCoeffMatrix('N', N, ...
                                    'Diffusion', Gamma / dx, ...
                                    'Convection', rho * U, ...
                                    'PhiLeft', phi_A, ...
                                    'PhiRight', phi_B);
    % Solve
    phi_numerical = A \ B;
    phi_FVM_Num_full = [phi_A; phi_numerical; phi_B];

    % % Analytical solution
    % numerator = exp((rho * U * x_analytical) / Gamma) - 1;
    % denominator = exp((rho * U * L) / Gamma) - 1;
    % phi_analytical = phi_A + (phi_B - phi_A) * (numerator / denominator);

     % Analytical solution
    phi_analytical = anaylitical_fun(x_analytical);

    % Plot
    figure;
    hold on;
    plot(x_analytical, phi_FVM_Num_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'FVM Numerical Solution');
    plot(x_analytical, phi_analytical, 'b--', 'LineWidth', 2, 'DisplayName', 'Analytical Solution');
    xlabel('Distance (m)', 'FontSize', 14);
    ylabel('\phi', 'FontSize', 14);
    title(sprintf('Numerical vs Analytical Solution (U = %.2f, N = %d)', U, N), 'FontSize', 14);
    legend('Location', 'best');
    grid on;
    hold off;

    % Print table
    Difference = abs(phi_FVM_Num_full - phi_analytical');
    Percentage_Error = (Difference ./ phi_analytical') * 100;

    T = table((0:N+1)', x_analytical', phi_FVM_Num_full, phi_analytical', Difference, Percentage_Error, ...
        'VariableNames', {'Node', 'Distance_m', 'FVM_Solution', 'Analytical_Solution', 'Difference', 'Percentage_Error'});
    fprintf(repmat('-', 1, 100));
    fprintf("\nThe No.Of Nodes: %d, Velocity: %.4f m/s\n", N, U);
    fprintf(repmat('-', 1, 20));
    fprintf('\n%-6s %-12s %-20s %-20s %-15s %-15s\n', 'Node', 'Distance(m)', 'FVM Solution', 'Analytical Solution', 'Difference', 'Percentage Error');
    fprintf(repmat('-', 1, 100));
    fprintf('\n');

    for k = 1:height(T)
        fprintf('%-4d %-12.4f %-20.6f %-20.4f %-15.3f %-15.2f\n', ...
                T.Node(k), T.Distance_m(k), T.FVM_Solution(k), ...
                T.Analytical_Solution(k), T.Difference(k), T.Percentage_Error(k));
    end
    disp('------ End of the Table ------------')
end


% Loop through each combination of U and N
for i = 1:length(U_values)
    for j = 1:length(N_values)
        U = U_values(i);  % Current velocity
        N = N_values(j);  % Current number of nodes
        SolveCDsConvectionDiffusionCase('U', U, 'N', N, 'phi_A', 1, 'phi_B', 0);
    end
end

