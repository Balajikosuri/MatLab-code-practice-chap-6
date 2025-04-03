function [A, B] = TriDiagonalCoeffMatrix(varargin)
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
    a_W = D + F/2;  % West coefficient
    a_E = D - F/2;  % East coefficient

    % Loop over each node in the grid
    for i = 1:N
        if i == 1  % Left Boundary Node (P = 1)
            A(i, i) = 3*D + (F/2);  % Modified coefficient for boundary
            A(i, i+1) = -a_E;  % East neighbor
            B(i) = (2*D + F) * phi_A;  % Apply left boundary flux

        elseif i == N  % Right Boundary Node (P = N)
            A(i, i) = 3*D - F/2;  % Modified coefficient for boundary
            A(i, i-1) = -a_W;  % West neighbor
            B(i) = (2*D - F) * phi_B;  % Apply right boundary flux

        else  % Internal Nodes (2 ≤ P ≤ N-1)
            A(i, i) = a_W + a_E - S_p(i);  % Central coefficient with source term
            A(i, i-1) = -a_W;  % West neighbor
            A(i, i+1) = -a_E;  % East neighbor
            B(i) = S_u(i);  % Source term contribution
        end
    end
end











% function [A, B] = TriDiagonlCoeffMatrix(N, a_W, a_E, T_A, T_B,S_p,S_u,q_R)
% 
%     %fron here start to
%     for i = 1:N
%         if i == 1  % Left Boundary Node (P = 1)
%             % ap = aw+aE-sp aW =(D-F/2) aE = (D+F/2)
%             % A(i, i) = (0 + (D-F/2)+ (2*D+F));
%             A(i, i) = 3*D-(F/2);
%             A(i, i+1) = -a_E;
%             B(i) = (2*D+F)*phi_A;  % Apply left boundary contribution 
%         elseif i == N  % Right Boundary Node (P = N)
%             % ap = aw+aE-sp
%             % A(i, i) = ((D-F/2)+0 +(2*D-F)); i.e
%             A(i,i) = 3*D-3*(F/2);
%             A(i, i-1) = -a_W;
%             B(i) = (2*D-F)*phi_A;  % Apply right boundary contribution
%         else  % Internal Nodes (2 ≤ P ≤ N-1)
%             A(i, i) = (a_W + a_E);
%             A(i, i-1) = -a_W;
%             A(i, i+1) = -a_E;
%             B(i) = 0;  % No source term in internal nodes
%         end
%     end
%     % here code every this fine but write the name value function 
% end
