function [phi] = QuickSolvePhiByGS(varargin)
    % QUICK-based 1D convection-diffusion solver using Gauss-Seidel
    % Uses name-value pairs
    % Derives dx = L/N, D = Gamma/dx, F = rho*U internally

    % Parse inputs
    p = inputParser;
    addParameter(p, 'N', []);
    addParameter(p, 'L', []);
    addParameter(p, 'Gamma', []);
    addParameter(p, 'rho', []);
    addParameter(p, 'U', []);
    addParameter(p, 'phi_A', []);
    addParameter(p, 'phi_B', []);
    addParameter(p, 'tol', 1e-6);  % Default tolerance
    parse(p, varargin{:});

    % Assign parsed values
    N = p.Results.N;
    L = p.Results.L;
    Gamma = p.Results.Gamma;
    rho = p.Results.rho;
    U = p.Results.U;
    phi_A = p.Results.phi_A;
    phi_B = p.Results.phi_B;
    tol = p.Results.tol;

    % Derived values
    dx = L / N;
    D = Gamma / dx;
    F = rho * U;

    % Initialization
    phi = zeros(N, 1);
    error = 1;

    while error > tol 
        phi_old = phi;

        for i = 1:N
            if i == 1
                aWW_eff = 0;
                aW_eff  = 0;
                aE_eff  = D + D/3 - 3*F/8;
                SP_eff  = 8*D/3 + 2*F/8 + F;
                SU_eff  = SP_eff * phi_A;
                aP      = aWW_eff + aW_eff + aE_eff + SP_eff;
                phi(i)  = (aE_eff * phi(i + 1) + SU_eff) / aP;

            elseif i == 2
                aWW_eff = 0;
                aW_eff  = D + F;
                aE_eff  = D - 3*F/8;
                SP_eff  = -F/4;
                SU_eff  = SP_eff * phi_A;
                aP      = aWW_eff + aW_eff + aE_eff + SP_eff;
                phi(i)  = (aW_eff * phi(i - 1) + aE_eff * phi(i + 1) + SU_eff) / aP;

            elseif i == N
                aWW_eff = -F/8;
                aW_eff  = 4*D/3 + 6*F/8;
                aE_eff  = 0;
                SP_eff  = 8*D/3 - F;
                SU_eff  = SP_eff * phi_B;
                aP      = aWW_eff + aW_eff + aE_eff + SP_eff;
                phi(i)  = (aWW_eff * phi(i - 2) + aW_eff * phi(i - 1) + SU_eff) / aP;

            else
                aWW_eff = -F/8;
                aW_eff  = D + 7*F/8;
                aE_eff  = D - 3*F/8;
                SP_eff  = 0;
                SU_eff  = 0;
                aP      = aWW_eff + aW_eff + aE_eff + SP_eff;
                phi(i)  = (aWW_eff * phi(i - 2) + aW_eff * phi(i - 1) + ...
                           aE_eff * phi(i + 1) + SU_eff) / aP;
            end
        end

        error = max(abs(phi - phi_old));
    end
end
