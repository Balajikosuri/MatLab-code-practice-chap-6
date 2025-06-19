% 2 D TDMA %
%%%% Variables definition %%%%%
clear;close all; clc;



nx = 30; % number of points in x;
ny = 3; % number of points in y;
L=0.2;
B=0.04;
gamma = 0.000217;
phi_top = 0;
phi_bottom = 40;


%%%%%
%Considering equal dx,dy and cell area%%
dx =L/nx;
dy =B/ny;
DifCof  = gamma /dy ;

%% initalize the nbhd all coeff 
aw=zeros(nx,ny)+(1/2);
ae=zeros(nx,ny)+(-1/2);
an=zeros(nx,ny)+(DifCof);
as=zeros(nx,ny)+(DifCof);
ap=zeros(nx,ny);
Sp=zeros(nx,ny);
Su=zeros(nx,ny);

Phi=zeros(nx,ny);

%   Discretized equation : aPTP = aWTW + aETE + aSTS + aNTN
%   Interior points : ap=aw+ae+as+an

for i = 2:nx-1
    for j = 2:nx
        ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j);
    end
end

%Boundaries: ap=aw+ae+as+an+Sp

for i = 1:nx
    if i==1 %Top boundary
        for j = 1:ny
            if j==1
                aw(i,j) = 0;
                an(i,j) = 0;
                Su(i,j)=2*DifCof*phi_top(1);
                Sp(i,j)=2*DifCof - 1/2;
                
            elseif j==ny
                % fprintf('balaji Line 58 (i,j) - (%f,%f)\n',i,j)
                ae(i,j)=0;
                an(i,j)=0;
                Su(i,j)=2*DifCof*phi_top;
                Sp(i,j)=2*DifCof+1/2;
                
            else 
                % fprintf('balaji Line 64 (i,j) - (%f,%f)\n',i,j)
                an(i,j)=0;
                Su(i,j)=2*DifCof*phi_bottom;
                Sp(i,j)=2*DifCof;
            end
        end
    elseif i==nx % bottom boundary
        for j = 1:ny
            if j==1
                aw(i,j)=0; as(i,j)=0;
                Su(i,j)=2*DifCof*phi_bottom;
                Sp(i,j)=-1/2;
            elseif j==ny 
                as(i,j)=0; ae(i,j)=0;
                Su(i,j)=2*DifCof*phi_bottom;
                Sp(i,j)=2*DifCof+1/2;    
            else 
                as(i,j)=0;
                Su(i,j)=2*DifCof*phi_bottom; 
                Sp(i,j)=2*DifCof;
            end
        end
    else 
        for j = 1:ny 
            if j==1
                % fprintf('balaji Line 90 (i,j)(%f,%f)\n',i,j)
                aw(i,j)=0;
                Sp(i,j)=-1/2;
                
            elseif j==ny  
                ae(i,j)=0;
                Sp(i,j)=1/2;
            end
        end
    end
    for j=1:ny
        ap(i,j)=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+Sp(i,j);
    end
end




function [T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi)
    T = Phi;
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

[T, run, max_res] = solve_tdma_rowwise(nx, ny, aw, ae, an, as, ap, Su, Phi);

Phi = T;
disp(Phi)
% Framing Temp Full 
% %%
% % Plot results at selected x locations
% x_sample = [0.2, 0.4, 0.6, 1.0];
% colors = 'rbgk';
% figure; hold on
% for k = 1:length(x_sample)
%     idx = round(x_sample(k) / dx) + 1;
%     disp(idx)
%     disp(Phi(:, idx))
%     plot(Phi(:, idx), y, colors(k), 'DisplayName', ['x = ' num2str(x_sample(k))])
% end
% xlabel('\phi'); ylabel('y'); grid on
% title('Convection-Diffusion: \phi(y) at Selected x')
% legend show
% %%


disp('T_interior')
disp(flip(Phi))
phiTop = phi_top;     % Top row fixed value
phiBottom = phi_bottom;  % Bottom row fixed value
T_full = frame_phi_full(Phi, phiTop, phiBottom);
disp(T_full);

% Framing X Y Mesh
[X, Y] = meshgrid(1:(nx+2), 1:(ny+2));

disp('Phi full')
disp(T_full)
% Plot Graph
contourf(X,Y,flipud(T_full),'ShowText','on')
colorbar; title('Temperature distribution (^o C)','FontSize',22);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);


% % Labels 
% xlabel('Bottom boundary Insulated:dT_B/dy = 0');
% ylabel('steady heat flux of 500x10^3 W/m^2');
% title('2D Plate - West: heatlux ; South & East: Insulated');
% subtitle('Top Boundary = 100 ^o C');

% Right-side label
% text(max(xlim), mean(ylim), '   Right boundary Insulated:dT_R/dx = 0', 'Rotation', 90, 'HorizontalAlignment', 'left');

% % printing the status of iteration 
% if max_res<(error)
%     fprintf('Solution converged');
% else
%     fprintf('Max. iterations reached');
% end


function Phi_full = frame_phi_full(Phi_interior, phiTop, phiBottom)
% Frames the interior Phi matrix with boundary conditions:
% - Top boundary: fixed value phiTop (scalar)
% - Bottom boundary: fixed value phiBottom (scalar)
% - Left boundary: copy of first interior column (succeeding)
% - Right boundary: copy of last interior column (preceding)
%
% Inputs:
%   Phi_interior : ny x nx matrix of interior values
%   phiTop       : scalar, value at top boundary row
%   phiBottom    : scalar, value at bottom boundary row
%
% Output:
%   Phi_full     : (ny+2) x (nx+2) matrix including boundaries

    [ny, nx] = size(Phi_interior);

    % Initialize full matrix
    Phi_full = zeros(ny + 2, nx + 2);

    % Fill interior
    Phi_full(2:end-1, 2:end-1) = Phi_interior;

    % Top boundary (first row) = phiTop scalar repeated
    Phi_full(1, :) = phiTop;

    % Bottom boundary (last row) = phiBottom scalar repeated
    Phi_full(end, :) = phiBottom;

    % Left boundary (first column excluding corners)
    Phi_full(2:end-1, 1) = Phi_interior(:, 1);

    % Right boundary (last column excluding corners)
    Phi_full(2:end-1, end) = Phi_interior(:, end);

    % Corners remain zero (or can be set as needed)
end

