nx=3; % number of points in x;
ny=3; % number of points in y;
T_initial=zeros(nx,ny);
Su = ...
  [0 0.05208 0;
   0 0 0;
   0.05208 0.05208 0.05208];

ae = ...
  [-0.5 -0.5 0;
   -0.5 -0.5 0;
   -0.5 -0.5 0];

an = ...
  [0 0 0;
   0.000651 0.000651 0.000651;
   0.000651 0.000651 0.000651];

ap = ...
  [-0.99804699999999991 0.0019529999999999999 1.0019529999999999;
   -0.998698 0.001302 1.001302;
   -0.999349 0.0019529999999999999 1.0019529999999999];

as = ...
  [0.000651 0.000651 0.000651;
   0.000651 0.000651 0.000651;
   0 0 0];

aw = ...
  [0 0.5 0.5;
   0 0.5 0.5;
   0 0.5 0.5];




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

[T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, T_initial);

% Display
disp('Final T (after convergence):');
disp(T);
