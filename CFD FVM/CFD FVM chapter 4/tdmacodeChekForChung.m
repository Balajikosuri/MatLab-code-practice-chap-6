clc; clear;

nx = 3; ny = 3; % Grid size
T_initial = zeros(ny, nx); % ny rows, nx columns

% Define coefficients (make sure they match (ny x nx) format)
Sp = [-0.498698 0.001302 0.501302;
      -0.5      0        0.5;
      -0.5      0.001302 0.501302];

Su = [0        0.05208 0;
      0        0       0;
      0.05208  0.05208 0.05208];

ae = [-0.5 -0.5 0;
      -0.5 -0.5 0;
      -0.5 -0.5 0];

an = [0        0        0;
      0.000651 0.000651 0.000651;
      0.000651 0.000651 0.000651];

ap = [-0.998047 0.001953 1.001953;
      -0.998698 0.001302 1.001302;
      -0.999349 0.001953 1.001953];

as = [0.000651 0.000651 0.000651;
      0.000651 0.000651 0.000651;
      0        0        0];

aw = [0 0.5 0.5;
      0 0.5 0.5;
      0 0.5 0.5];

% Print all coefficients row by row
for i = 1:ny
    fprintf('--- Row %d ---\n', i);
    for j = 1:nx
        fprintf('aw(%d,%d) = %.6f, ae = %.6f, an = %.6f, as = %.6f, ap = %.6f, Sp = %.6f, Su = %.6f\n', ...
            i, j, aw(i,j), ae(i,j), an(i,j), as(i,j), ap(i,j), Sp(i,j), Su(i,j));
    end
end

% Solve
[T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, T_initial);

% Display
disp('Final T (after convergence):');
disp(T);

% --- Row-wise TDMA Solver ---
function [T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, T_initial)
    T = T_initial;
    run = 1;
    max_res = 10;
    tol = 1e-3;

    while max_res > tol
        T_old = T;

        % Loop over each row (i.e., fixed y, solve in x-direction)
        for i = 1:ny
            alpha = zeros(1, nx);
            beta = zeros(1, nx);
            D = zeros(1, nx);
            C = zeros(1, nx);
            A = zeros(1, nx);
            Cdash = zeros(1, nx);

            for j = 1:nx
                alpha(j) = ae(i, j);
                beta(j) = aw(i, j);
                D(j) = ap(i, j);
                C(j) = Su(i, j);

                % Add vertical neighbors' contributions
                if i > 1
                    C(j) = C(j) + as(i, j) * T(i - 1, j);
                end
                if i < ny
                    C(j) = C(j) + an(i, j) * T(i + 1, j);
                end
            end

            % Forward sweep
            for j = 1:nx
                if j == 1
                    A(j) = alpha(j) / D(j);
                    Cdash(j) = C(j) / D(j);
                else
                    denom = D(j) - beta(j) * A(j - 1);
                    A(j) = alpha(j) / denom;
                    Cdash(j) = (C(j) + beta(j) * Cdash(j - 1)) / denom;
                end
            end

            % Backward substitution
            for j = nx:-1:1
                if j == nx
                    T(i, j) = Cdash(j);
                else
                    T(i, j) = A(j) * T(i, j + 1) + Cdash(j);
                end
            end
        end

        run = run + 1;
        max_res = max(abs(T(:) - T_old(:)));
    end
end
