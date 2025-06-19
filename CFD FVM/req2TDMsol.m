clear;close all;clc;



% Sp = ...
%   [0 0 0 -20;
%    0 0 0 -20;
%    0 0 0 -20];
% 
% Su = ...
%   [500 500 500 2500;
%    0 0 0 2000;
%    0 0 0 2000];
% 
% an = ...
%   [10 10 10 0;
%    10 10 10 0;
%    10 10 10 0];
% 
% ap = ...
%   [20 30 30 40;
%    30 40 40 50;
%    20 30 30 40];
% 
% ae =[ 10	10	10	10;
% 10	10	10	10;
% 0	0	0	0  ; 
% ];
% 
% as = ...
%   [0 10 10 10;
%    0 10 10 10;
%    0 10 10 10];
% 
% aw = ...
%   [0 0 0 0;
%    10 10 10 10;
%    10 10 10 10];
% 
% Define 2x2 layout for each quantity using node map:
% [1 2
%  3 4]
% 
% Sp = [3, 1;
%       11, 9;];
% 
% Su = [310, 10;
%       910, 610];
% 
% ae = [0.5, 0;
%       0.5, 0];
% 
% aw = [0, 1.5;
%       0, 1.5];
% 
% an = [0, 0;
%       -1, -1];
% 
% as = [3, 3;
%       0, 0];
% 
% ap = [6.5, 5.5;
%       10.5, 9.5];


function [aw, ae, an, as, ap, Sp, Su] = apply_boundary_conditions(nx, ny, aw, ae, an, as, ap, Sp, Su)
    % Applies boundary conditions to the coefficient arrays.

    % Interior points : ap=aw+ae+as+an
    for i = 2:nx - 1
        for j = 2:ny - 1 % Corrected inner loop limit for interior
            ap(i, j) = aw(i, j) + ae(i, j) + as(i, j) + an(i, j);
        end
    end

    % Boundaries: ap=aw+ae+as+an-Sp
    for i = 1:nx
        if i == 1 % left boundary
            for j = 1:ny
                if j == 1
                    aw(i, j) = 0;
                    as(i, j) = 0;
                    Su(i, j) = 500;
                    Sp(i, j) = 0;

                elseif j == ny
                    aw(i, j) = 0;
                    an(i, j) = 0;
                    Su(i, j) = 2500;
                    Sp(i, j) = -20;

                else
                    aw(i, j) = 0;
                    Su(i, j) = 500;
                    Sp(i, j) = 0;
                end
            end
        elseif i == nx % right boundary
            for j = 1:ny
                if j == 1
                    ae(i, j) = 0;
                    as(i, j) = 0;
                    Su(i, j) = 0;
                    Sp(i, j) = 0;

                elseif j == ny
                    an(i, j) = 0;
                    ae(i, j) = 0;
                    Su(i, j) = 2000;
                    Sp(i, j) = -20;

                else
                    Su(i, j) = 0;
                    ae(i, j) = 0;
                    Sp(i, j) = 0;
                end
            end
        else % interior columns (not left or right boundary)
            for j = 1:ny
                if j == 1 % bottom boundary in interior columns
                    as(i, j) = 0;
                    Su(i, j) = 0;
                    Sp(i, j) = 0;

                elseif j == ny % top boundary in interior columns
                    an(i, j) = 0;
                    Su(i, j) = 2000;
                    Sp(i, j) = -20;
                end
                % Note: Interior points in interior columns already have their ap set
                % This section only adjusts boundary points in interior columns
            end
        end
        % Recalculate ap for all points in the current column after boundary adjustments
        for j = 1:ny
            ap(i, j) = aw(i, j) + ae(i, j) + as(i, j) + an(i, j) - Sp(i, j);
        end
    end
end

nx=3; % number of points in x;
ny=3; % number of points in y;
T_initial=zeros(nx,ny);
Sp = ...
  [-0.498698 0.001302 0.501302;
   -0.5 0 0.5;
   -0.5 0.001302 0.501302];

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


function [T, run, max_res] = solve_tdma(nx, ny, aw, ae, an, as, ap, Su, T_initial)
    % Solves the system of equations using TDMA.
    T = T_initial; % Start with the initial temperature guess

    alpha = zeros(nx, ny);
    beta = zeros(nx, ny);
    D = zeros(nx, ny);
    C = zeros(nx, ny);

    A = zeros(nx, ny);
    Cdash = zeros(nx, ny);

    run = 1;
    error = 0.001; % set desired convergence value
    max_res = 10; % Initialize max_res high to enter the loop
    while max_res > error 
        t = T; % Store current T for residual calculation

        for i = 1:nx
            for j = 1:ny
                alpha(i, j) = an(i, j);
                beta(i, j) = as(i, j);
                D(i, j) = ap(i, j);

                % Calculate C based on neighboring temperatures and source terms
                if i == 1 % left boundary column
                     if nx > 1 % Ensure there is a point to the east
                         C(i, j) = (ae(i, j) * T(i + 1, j)) + Su(i, j);
                     else % Handle case of a 1D problem in x if nx=1
                         C(i, j) = Su(i, j);
                     end
                elseif i == nx % right boundary column
                    if nx > 1 % Ensure there is a point to the west
                         C(i, j) = (aw(i, j) * T(i - 1, j)) + Su(i, j);
                     else % Handle case of a 1D problem in x if nx=1
                         C(i, j) = Su(i, j);
                     end
                else % interior columns
                    C(i, j) = (aw(i, j) * T(i - 1, j)) + (ae(i, j) * T(i + 1, j)) + Su(i, j);
                end
            end

            % Forward elimination (across rows, for a fixed i)
            for j = 1:ny
                if j == 1
                    A(i, j) = (alpha(i, j) / (D(i, j)));
                    Cdash(i, j) = (C(i, j) / (D(i, j)));
                else
                    A(i, j) = (alpha(i, j) / (D(i, j) - (beta(i, j) * A(i, j - 1))));
                    Cdash(i, j) = (((beta(i, j) * Cdash(i, j - 1)) + C(i, j)) / (D(i, j) - (beta(i, j) * A(i, j - 1))));
                end
            end

            % Backward substitution (across rows, for a fixed i)
            for j = ny:-1:1
                if j == ny
                    T(i, j) = Cdash(i, j);
                else
                    T(i, j) = (A(i, j) * T(i, j + 1)) + Cdash(i, j);
                end
            end
        end

        run = run + 1;
        res = abs(t - T);
        max_res = max(res(:)); % Find max residual in the entire matrix

        % fprintf('Iteration.... %i\n', run); % Uncomment to see iteration progress
    end
end


[T, run, max_res] = solve_tdma(nx, ny, aw, ae, an, as, ap, Su, T_initial);



disp(flip(T'))