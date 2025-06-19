clc; clear; close all;
T_right = 0; % at t>0 

% Initial temperature at t = 0
T_old = [200; 200; 200; 200; 200];

% Coefficient matrix A
A = [2125   -125     0     0       0;
     -125   2250   -125    0       0;
       0   -125     2250  -125     0;
       0     0     -125    2250  -125;
       0     0      0     -125  2375];

% Time step info (adjust if needed)
dt = 2;
t_total = [40, 80, 120];  % output times
steps_needed = t_total / dt;

% Store results
T_all = zeros(5, length(steps_needed));

% Time loop
for step = 1:max(steps_needed)
    % Construct RHS using your custom formula
    b = ones(5,1);
    b(1) = 1875*T_old(1) + 125*T_old(2);
    b(2) = 125*T_old(1) + 1750*T_old(2) + 125*T_old(3);
    b(3) = 125*T_old(2) + 1750*T_old(3) + 125*T_old(4);
    b(4) = 125*T_old(3) + 1750*T_old(4) + 125*T_old(5);
    b(5) = 125*T_old(4) + 1625*T_old(5)+250*T_right;
    
    % Solve linear system A * T_new = b
    T_new = A \ b;

    % Store and display results at desired time steps
    match = find(step == steps_needed);
    if ~isempty(match)
        T_all(:, match) = T_new;
        fprintf('Time = %d s\n', step * dt);
        disp(T_new);
    end

    % Update for next time step
    T_old = T_new;
end