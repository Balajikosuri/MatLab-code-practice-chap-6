function [A, B] = TriDiagonlCoeffMatrix(N, a_W, a_E, T_A, T_B,S_p,S_u,q_R)
    % Initialize Coefficient Matrix A and RHS Vector B
    A = zeros(N, N);
    B = zeros(N, 1);
    if nargin < 6  % nargin gives the number of input arguments provided
        S_p= 0; % Set S_p to zero if not provided
    end
    if nargin < 7  % nargin gives the number of input arguments provided
        S_u = 0; % Set S_u to zero if not provided
    end
  
    if nargin < 8  % nargin gives the number of input arguments provided
        q_R= NaN; % Set q to NA if not provided
    end
    % Construct Coefficient Matrix A and RHS Vector B using For Loop
    for i = 1:N
        if i == 1  % Left Boundary Node (P = 1)
            A(i, i) = (2 * a_W + a_E + S_u);
            A(i, i+1) = -a_E;
            B(i) = S_p + 2 * a_W * T_A;  % Apply left boundary contribution
        elseif i == N  % Right Boundary Node (P = N)
            if q_R == 0
                A(i, i) = (a_W + S_u);
                A(i, i-1) = -a_W;
                B(i) = S_p;
            else
                A(i, i) = (a_W + 2 * a_E + S_u);
                A(i, i-1) = -a_W;
                B(i) = S_p + 2 * a_E * T_B;  % Apply right boundary contribution
            end
        else  % Internal Nodes (2 ≤ P ≤ N-1)
            A(i, i) = (a_W + a_E+S_u);
            A(i, i-1) = -a_W;
            A(i, i+1) = -a_E;
            B(i) = S_p;  % No source term in internal nodes
        end
    end
end
