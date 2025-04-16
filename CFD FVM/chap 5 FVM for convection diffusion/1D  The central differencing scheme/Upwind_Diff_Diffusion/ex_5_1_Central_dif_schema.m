clc; clear; close all;

% Given Parameters
L = 1.0;        % Length of the domain
Gamma = 0.1;    % Diffusion coefficient (Γ)
rho = 1.0;      % Density

phi_A = 1; % Left boundary condition
phi_B = 0; % Right boundary condition

% Define different values of U (velocity) and N (grid points)
U_values = [0.1, 2.5];  % Cases: Low and high velocity
N_values = [5,];     % Cases: Coarse and fine grid

% Loop through each combination of U and N
for i = 1:length(U_values)
    for j = 1:length(N_values)
        U = U_values(i);  % Current velocity
        N = N_values(j);  % Current number of nodes
        dx = L / N; % Grid spacing
        % Compute x_analytical values (cell-centered positions)
        P = 1:N; % Node indices
        x_P = (P - 0.5) * dx;  % Apply x_P = x_0 + (P - 1/2) * dx

        % Include boundary points
        x_analytical = [0, x_P, L];
        % Generate coefficient matrix and source term
        [A, B] = TriDiagonalCoeffMatrix('N', N, 'Diffusion', Gamma / dx, ...
                                        'Convection', rho * U, 'PhiLeft', phi_A, 'PhiRight', phi_B);
        
        % Solve for phi
        phi_numerical = A \ B;
        
        % Compute analytical solution
        numerator = exp((rho*U*x_analytical)/Gamma)-1;
        denominator = exp((rho*U*L)/Gamma)-1;
        phi_analytical = phi_A+(phi_B-phi_A).*(numerator/denominator);
        % disp(phi_analytical)
        % Display computed x_analytical values for verification
        phi_FVM_Num_full = [phi_A;phi_numerical;phi_B];
        % Plot for this case
        figure;
        hold on;
        plot(x_analytical, phi_FVM_Num_full, 'rs-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'FVM Numerical Solution');
        plot(x_analytical, phi_analytical, 'b--', 'LineWidth', 2, 'DisplayName', 'Analytical Solution');
        
        xlabel('Distance (m)', 'FontSize', 14);
        ylabel('\phi', 'FontSize', 14);
        title(sprintf('Numerical vs Analytical Solution (U = %.1f, N = %d)', U, N), 'FontSize', 14);
        legend('Location', 'best');
        grid on;
        hold off;
        % Printing Table 
        % Compute differences and percentage errors
        Difference = abs(phi_FVM_Num_full - phi_analytical');
        Percentage_Error = (Difference ./ phi_analytical') * 100;
        
        % Create a table
        T = table((0:N+1)', x_analytical', phi_FVM_Num_full, phi_analytical', Difference, Percentage_Error, ...
            'VariableNames', {'Node', 'Distance_m', 'FVM_Solution', 'Analytical_Solution', 'Difference', 'Percentage_Error'});
        
        % Display table
        % disp(T);
        
        % Print table with formatted output
        fprintf("\nThe No.Of Nodes: %d velocity: %.4f m/s  \n",N,U)
        fprintf('\n%-6s %-12s %-20s %-20s %-15s %-15s\n', 'Node', 'Distance(m)', 'FVM Solution', 'Analytical Solution', 'Difference', 'Percentage Error');
        fprintf(repmat('-', 1, 100));
        fprintf('\n');
        
        for k = 1:height(T)
            fprintf('%-4d %-12.4f %-20.4f %-20.4f %-15.3f %-15.2f\n', T.Node(k), T.Distance_m(k), ...
                T.FVM_Solution(k), T.Analytical_Solution(k), T.Difference(k), T.Percentage_Error(k));
        end
        disp('------ End of the Table ------------')
        % Tablen End;
       
    end
end
