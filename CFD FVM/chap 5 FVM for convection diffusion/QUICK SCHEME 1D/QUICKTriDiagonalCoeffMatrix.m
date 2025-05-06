function [A, B] = UDSTriDiagonalCoeffMatrix(varargin)
    % TriDiagonalCoeffMatrix: Generates the coefficient matrix A and source vector B
    % for a 1D convection-diffusion problem using the central differencing scheme.
    %
    % Name-Value Pair Arguments:
    % 'N'         - Number of grid points (default: 5)
    % 'Diffusion' (D)  - Diffusion coefficient per control volume (default: 0.1)
    % 'Convection' (F) - Convective flux per control volume (default: 0.1)
    % 'PhiLeft'  (phi_A) - Left boundary condition (default: 1)
    % 'PhiRight' (phi_B) - Right boundary condition (default: 0)
    % 'SourceP'  (S_p) - Source coefficient (default: zeros(N,1))
    % 'SourceU'  (S_u) - Source term (default: zeros(N,1))
    %
    % OUTPUTS:
    % A - Tridiagonal coefficient matrix (NxN)
    % B - Source term vector (Nx1)

    % Parse Name-Value Arguments
    p = inputParser;
    addParameter(p, 'N', 5, @(x) isnumeric(x) && x > 0); % Default N = 5
    addParameter(p, 'Diffusion', 0.1, @isnumeric);
    addParameter(p, 'Convection', 0.1, @isnumeric);
    addParameter(p, 'PhiLeft', 1, @isnumeric);
    addParameter(p, 'PhiRight', 0, @isnumeric);
    
    % Parse inputs
    parse(p, varargin{:});
    
    % Assign parsed values to variables
    N = p.Results.N;
    D = p.Results.Diffusion;
    F = p.Results.Convection;
    phi_A = p.Results.PhiLeft;
    phi_B = p.Results.PhiRight;
    % Default source terms (zero by default)
    % S_p = zeros(N, 1);
    % S_u = zeros(N, 1);

    % Initialize matrices
    A = zeros(N, N);
    B = zeros(N, 1);
    D_A = D;
    D_B = D;

    % % Define coefficients
    % a_WW = -F/8; % mirror node 
    % a_W = (D + D_A/3 + (6*F)/8);  % West coefficient
    % a_E = ((5*D)/3 - (3*F)/8);  % East coefficient
    % %internal node coeff
    % a_W_int = D + (6*F)/8 + F/8;
    % a_E_int = D - (3*F)/8;
    % a_WW_int = -F/8 ;
    % a_p_int = a_W_int + a_E_int + a_WW_int;
 

    % Loop over each node in the grid
    for i = 1:N
        if i == 1  % Left Boundary Node (P = 1)
            a_WW = 0;
            a_W= 0;
            a_E = D + D_A/3 - 3*F/8;
            S_P = 8*D_A/3 + 2*F/8 + F;
            S_U = S_P*phi_A;
            A(i, i) = a_WW + a_W + a_E + S_P;  % aWW+aW+aE+Sp coefficient for boundary///
            A(i, i+1) = - a_E;  % East neighbor
            B(i) = S_U;  % Apply left boundary flux
        elseif i == 2
            a_WW = 0;
            a_W= D+F;
            a_E = D - 3*F/8;
            S_P = -F/4 ;
            S_U = S_P*phi_A;
            A(i, i) = a_WW + a_W + a_E + S_P;  % aWW+aW+aE+Sp coefficient for boundary///
            A(i, i+1) = - a_E;  % East neighbor
            A(i, i-1) = - a_W ; % west node coeff
            B(i) = S_U;  % Apply left boundary flux
        elseif i == N  % Right Boundary Node (P = N)
            a_WW = -F/8;
            a_W= 4*D/3 + 6*F/8;
            a_E = 0;
            S_P = 8*D/3 - F ;
            S_U = S_P*phi_B;
            A(i, i) = a_WW + a_W + a_E + S_P;  % aWW+aW+aE+Sp coefficient for boundary///
            A(i, i-1) = - a_W;  % west neighbor
            A(i, i-2) = - a_WW ; % imaginary node coeff
            B(i) = S_U;  % Apply left boundary flux

        else  % Internal Nodes (3 ≤ P ≤ N-1)
            a_WW = -F/8;
            a_W= D + 7*F/8;
            a_E = D - 3*F/8;
            S_P = 0 ;
            S_U = 0;
            A(i, i) = a_WW + a_W + a_E + S_P;  % aWW+aW+aE+Sp coefficient for boundary///
            A(i, i+1) = - a_E;  % East neighbor
            A(i, i-1) = - a_W;  % west neighbor
            A(i, i-2) = - a_WW ; % imaginary node coeff

            B(i) = S_U;  % Apply left boundary flux

        end
    end
end