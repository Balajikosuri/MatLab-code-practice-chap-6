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
    S_p = zeros(N, 1);
    S_u = zeros(N, 1);

    % Initialize matrices
    A = zeros(N, N);
    B = zeros(N, 1);

    % Define coefficients
    a_W = D + F;  % West coefficient
    a_E = D;  % East coefficient

    % Loop over each node in the grid
    for i = 1:N
        if i == 1  % Left Boundary Node (P = 1)
            A(i, i) = 3*D+F;  % Modified coefficient for boundary
            A(i, i+1) = -a_E;  % East neighbor
            B(i) = (2*D + F) * phi_A;  % Apply left boundary flux

        elseif i == N  % Right Boundary Node (P = N)
            A(i, i) = 3*D+F;  % Modified coefficient for boundary
            A(i, i-1) = -a_W;  % West neighbor
            B(i) = 2*D*phi_B;  % Apply right boundary flux

        else  % Internal Nodes (2 ≤ P ≤ N-1)
            A(i, i) = a_W + a_E ;  % Central coefficient with source term
            A(i, i-1) = -a_W;  % West neighbor
            A(i, i+1) = -a_E;  % East neighbor
            B(i) = 0;  % Source term contribution+
        end
    end
end