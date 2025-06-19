% clc;
% clear;
% 
% n = 5;
% 
% % Coefficients for the system
% a = [0, 0.55, 0.55, 0.55, 0.55];    % sub-diagonal (aW)
% b = [1.55, 1.0, 1.0, 1.0, 1.45];    % main diagonal (aP)
% c = [0.45, 0.45, 0.45, 0.45, 0];    % super-diagonal (aE)
% d = [1.1, 0, 0, 0, 0];              % RHS (Su), with φA = 1, φB = 0
% 
% % TDMA forward elimination
% for i = 2:n
%     m = a(i) / b(i-1);
%     b(i) = b(i) - m * c(i-1);
%     d(i) = d(i) - m * d(i-1);
% end
% 
% % TDMA backward substitution
% phi = zeros(n,1);
% phi(n) = d(n) / b(n);
% for i = n-1:-1:1
%     phi(i) = (d(i) - c(i) * phi(i+1)) / b(i);
% end
% 
% % Display results
% for i = 1:n
%     fprintf('φ%d = %.4f\n', i, phi(i));
% end
% 
% % %  
% 
% clc;
% clear;
% 
% % Number of nodes
% n = 5;
% 
% % Boundary values
% phi_A = 1.0;
% phi_B = 0.0;
% 
% % Coefficients from table
% a = [0, 2.5, 2.5, 2.5, 2.5];      % aW
% b = [3.5, 2.5, 2.5, 2.5, 3.5];    % aP = aW + aE - Sp
% c = [0, 0, 0, 0, 0];              % aE (all zero)
% d = [3.5 * phi_A, 0, 0, 0, 1.0 * phi_B];  % Su
% 
% % Thomas algorithm: forward elimination
% for i = 2:n
%     m = a(i) / b(i-1);
%     b(i) = b(i) - m * c(i-1);   % No effect as c = 0
%     d(i) = d(i) - m * d(i-1);
% end
% 
% % Back substitution
% phi = zeros(n, 1);
% phi(n) = d(n) / b(n);
% for i = n-1:-1:1
%     phi(i) = (d(i) - c(i) * phi(i+1)) / b(i); % c = 0 → simplified
% end
% 
% % Display results
% for i = 1:n
%     fprintf('phi(%d) = %.4f\n', i, phi(i));
% end
% 



%%% -----------------------------

clc;
clear;

% Number of nodes
n = 5;

% Boundary conditions
phi_A = 1.0;
phi_B = 0.0;

% % Coefficients from the table
% a = [0, 2.5, 2.5, 2.5, 2.5];         % aW
% b = [3.5, 2.5, 2.5, 2.5, 3.5];       % aP = aW + aE - Sp
% c = [0, 0, 0, 0, 0];                 % aE (all zero)
% d = [3.5 * phi_A, 0, 0, 0, 1.0 * phi_B];  % Su vector
a = [0, 0.55, 0.55, 0.55, 0.55];    % sub-diagonal (aW)
b = [1.55, 1.0, 1.0, 1.0, 1.45];    % main diagonal (aP)
c = [0.45, 0.45, 0.45, 0.45, 0];    % super-diagonal (aE)
d = [1.1, 0, 0, 0, 0];              % RHS (Su), with φA = 1, φB = 0

% Thomas algorithm - Forward elimination
for i = 2:n
    m = a(i) / b(i-1);
    b(i) = b(i) - m * c(i-1);       % c = 0 so this doesn't affect
    d(i) = d(i) - m * d(i-1);
end

% Backward substitution
phi = zeros(n, 1);
phi(n) = d(n) / b(n);
for i = n-1:-1:1
    phi(i) = (d(i) - c(i) * phi(i+1)) / b(i);  % c = 0 so simplified
end

% Display result
for i = 1:n
    fprintf('phi(%d) = %.4f\n', i, phi(i));
end

%%%%%%%%%%%%%%%%%

