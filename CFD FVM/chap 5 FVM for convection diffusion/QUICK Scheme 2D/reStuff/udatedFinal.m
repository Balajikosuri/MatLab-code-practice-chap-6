clear; clc; close all;
% Physical and boundary values
Lx = 3;     Ly = 3;     
Nx = 5;     Ny = 5;     
dx = Lx / Nx;   dy = Ly / Ny;

Phi_Left = 100; Phi_Right = 0;
Phi_Top = 0;    Phi_Bottom = 100;

rouh = 1;    
Gamma = 1;
a = 10;
b = 2;
u = 1;      
v = 4;
[phi, iter, error] = solveConvDiff2DimByCentralDS( ...
    'Lx', Lx, 'Ly', Ly, ...
    'Nx', Nx, 'Ny', Ny, ...
    'Phi_Left', Phi_Left, 'Phi_Right', Phi_Right, ...
    'Phi_Top', Phi_Top, 'Phi_Bottom', Phi_Bottom, ...
    'rouh', rouh, 'Gamma', Gamma, ...
    'a', a, 'b', b, ...
    'u', u, 'v', v, ...
    'tol', 1e-6);  

plotPhiResultsIn2Dim3Dim( ...
    'phi', phi, ...
    'dx', dx, 'dy', dy, ...
    'Nx', Nx, 'Ny', Ny, ...
    'Phi_Left', Phi_Left, 'Phi_Right', Phi_Right, ...
    'Phi_Top', Phi_Top, 'Phi_Bottom', Phi_Bottom, ...
    'error', error);


%% function for solving 2D Central Diff Scheme
function [phi, iterations, error] = solveConvDiff2DimByCentralDS(varargin)
% solveConvDiff2DimByCentralDS solves 2D convection-diffusion using FVM with Gauss-Seidel
%
% Usage:
%   [phi, iterations, error] = solveConvDiff2DimByCentralDS('Lx',1,'Ly',1,'Nx',20,'Ny',20,...)
%
% Key-Value Inputs:
%   'Lx', 'Ly'            - Length of domain (m)
%   'Nx', 'Ny'            - Number of control volumes
%   'Phi_Left', etc.      - Boundary values
%   'rouh', 'Gamma', 'a', 'b', 'u', 'v' - Fluid and source properties
%   'tol'                 - Convergence tolerance

    % Parse inputs
    p = inputParser;
    addParameter(p, 'Lx', 1);
    addParameter(p, 'Ly', 1);
    addParameter(p, 'Nx', 10);
    addParameter(p, 'Ny', 10);
    addParameter(p, 'Phi_Left', 0);
    addParameter(p, 'Phi_Right', 0);
    addParameter(p, 'Phi_Top', 0);
    addParameter(p, 'Phi_Bottom', 0);
    addParameter(p, 'rouh', 1);
    addParameter(p, 'Gamma', 1);
    addParameter(p, 'a', 0);
    addParameter(p, 'b', 0);
    addParameter(p, 'u', 0);
    addParameter(p, 'v', 0);
    addParameter(p, 'tol', 1e-6);
    parse(p, varargin{:});
    
    % Extract variables
    Lx = p.Results.Lx;
    Ly = p.Results.Ly;
    Nx = p.Results.Nx;
    Ny = p.Results.Ny;
    Phi_Left = p.Results.Phi_Left;
    Phi_Right = p.Results.Phi_Right;
    Phi_Top = p.Results.Phi_Top;
    Phi_Bottom = p.Results.Phi_Bottom;
    rouh = p.Results.rouh;
    Gamma = p.Results.Gamma;
    a = p.Results.a;
    b = p.Results.b;
    u = p.Results.u;
    v = p.Results.v;
    tol = p.Results.tol;

    % (Rest of your original function here, no changes needed)
    % ----------------------------------------------------------------------
    
    % Grid parameters
    dx = Lx / Nx;
    dy = Ly / Ny;
    
    % Coefficients for convection and diffusion
    Fe = rouh*u; Fw = rouh*u;
    Fn = rouh*v; Fs = rouh*v;

    D = Gamma / dx;
    De = D;
    Dw = D;
    Dn = D;
    Ds = D;
    S_u= 0; S_P = 0 ;

    % Initialize temperature field with boundary conditions
    phi = zeros(Ny + 2, Nx + 2);  % Including ghost cells
    phi(:, 1) = Phi_Left;          % Left boundary (x = 0)
    phi(:, end) = Phi_Right;       % Right boundary (x = Lx)
    phi(1, :) = Phi_Top;           % Top boundary (y = Ly)
    phi(end, :) = Phi_Bottom;      % Bottom boundary (y = 0)
    
    % Iterative solver (Gauss-Seidel)
    error = 1;
    iter = 0;
    % % Coefficients
    % aE_eff = De - (3/8)*Fe;
    % aW_eff = Dw + (6/8)*Fw;
    % aWW_eff = -(1/8)*Fw;
    % aN_eff = Dn - (3/8)*Fn;
    % aS_eff = Ds + (6/8)*Fs;
    % aSS_eff = -(1/8)*Fs;
    % % Modify coefficients near boundaries
    % aW_bdry_nxt = 2 * D_w + F_w / 2;
    % aE_bdry_nxt = 2 * D_e - F_e / 2;
    % aN_bdry_nxt = 2 * D_n - F_n / 2;
    % aS_bdry_nxt = 2 * D_s + F_s / 2;

    while error > tol
    phi_old = phi;

    % Update interior points
    for i = 2:Nx + 1
        for j = 2:Ny + 1
            
            % Initialize all coefficients and source terms
            aWW_eff = -Fw/8; aW_eff = Dw + 6*Fw/8 + Fe/8;
            aE_eff  = De - 3*Fe/8 ;
            aN_eff  = Dn - 3*Fn/8 ; 
            aSS_eff = -Fs/8; aS_eff = Ds + 6*Fs/8 + Fn/8; 
            S_P = 0;  S_u = 0;

            % %% X-direction: Left and Right boundaries
            % if i == 2  % Left boundary
            %     fprintf('line 142 i,j =>{%d,%d} \n',j,i)
            %     aW_eff = 0;
            %     S_P = b*dx*dy;
            % 
            % elseif i == 3 % nxet to boundary node genrally 2nd Node point (Left to Right) 
            %     aW_eff = Dw + 7*Fw/8 + Fe/8;
            %     aE_eff = De - 3*Fe/8;
            %     S_P = -Fw/4;
            %     S_u = (-Fw/4) * Phi_Left;
            % elseif i == Nx + 1  % Right boundary
            %     aW_eff = Dw + 6*Fw/8 + Fe/8;
            %     S_u  = S_P * Phi_Right;
            %     aE_eff = 0;
            % else  % Internal nodes
            %     aW_eff = D + 7*Fw/8 + Fe/8;
            %     aE_eff = Dw + 6*Fw/8 + Fe/8;
            % end
            % 
            % %% Y-direction: Top and Bottom boundaries
            % if j == 2  % Top boundary
            %     aN_eff = 0;
            % else
            %     aN_eff = Dw + 6*Fw/8 + Fe/8;
            % end
            % 
            % if j == Ny + 1  % Bottom boundary
            %     aS_eff = 0;
            % elseif i == 3
            %     aS_eff = Ds + 7*Fs/8 + Fn/8;
            % else
            %     aS_eff = Ds + 6*Fs/8 + Fn/8;
            % end
            % 
            %% Final coefficient assembly and update phi
            aP = aWW_eff + aW_eff + aE_eff + aN_eff + aS_eff + S_P;

            % Use different stencil depending on node location
            if i == 2 || j == Ny + 1
                phi(j, i) = (aE_eff * phi(j, i + 1) + aW_eff * phi(j, i - 1) ...
                           + aN_eff * phi(j - 1, i) + aS_eff * phi(j + 1, i) + S_u) / aP;
            else
                phi(j, i) = (aE_eff * phi(j, i + 1) + aW_eff * phi(j, i - 1) + aWW_eff * phi(j, i - 2) ...
                           + aN_eff * phi(j - 1, i) + aS_eff * phi(j + 1, i) + aSS_eff * phi(j + 2, i) + S_u) / aP;
            end
        end
    end

    % Compute error
    error = max(max(abs(phi - phi_old)));
    iter = iter + 1;
    end


    % Return the final solution and the number of iterations
    iterations = iter;

end

%% function for Ploting 2D & 3D 

function plotPhiResultsIn2Dim3Dim(varargin)
    % This function plots the solution of the 2D convection-diffusion problem.
    % It includes both 2D contour and 3D surface plots with boundary labels.
    % 
    % Inputs (key-value pairs):
    %   'phi'        - Full phi matrix (with ghost cells)
    %   'dx', 'dy'   - Cell dimensions
    %   'Nx', 'Ny'   - Number of cells in x and y directions
    %   'Phi_Left', 'Phi_Right', 'Phi_Top', 'Phi_Bottom' - Boundary phi values
    %   'error'      - Final error value

    % Parse Inputs
    p = inputParser;
    addParameter(p, 'phi', []);
    addParameter(p, 'dx', 0);
    addParameter(p, 'dy', 0);
    addParameter(p, 'Nx', 0);
    addParameter(p, 'Ny', 0);
    addParameter(p, 'Phi_Left', 0);
    addParameter(p, 'Phi_Right', 0);
    addParameter(p, 'Phi_Top', 0);
    addParameter(p, 'Phi_Bottom', 0);
    addParameter(p, 'error', 0);
    
    parse(p, varargin{:});
    
    phi = p.Results.phi;
    dx = p.Results.dx;
    dy = p.Results.dy;
    Nx = p.Results.Nx;
    Ny = p.Results.Ny;
    Phi_Left = p.Results.Phi_Left;
    Phi_Right = p.Results.Phi_Right;
    Phi_Top = p.Results.Phi_Top;
    Phi_Bottom = p.Results.Phi_Bottom;
    error = p.Results.error;

    %% 2D Contour Plot
    figure;
    contourf(flipud(phi), 12, 'LineColor', 'none'); 
    colorbar;
    colormap(jet); 
    title('2D Heat Conduction - Finite Volume Method 2D Upwind.DS');
    xlabel('x (m)');
    ylabel('y (m)');

    %% Wall Boundary Labels
    annotation('textbox', [0.43 0.97 0.13 0.03], ...
        'String', {sprintf('\\phi_{Top} = %g', Phi_Top)}, ...
        'EdgeColor','none', 'FontSize', 12, 'FontWeight', 'bold');
    
    annotation('textbox', [0.01 0.49 0.07 0.07], ...
        'String', {sprintf('\\phi_{Left} = %g', Phi_Left)}, ...
        'EdgeColor','none', 'FontSize', 12, 'FontWeight', 'bold');
    
    annotation('textbox', [0.44 0.01 0.14 0.04], ...
        'String', {sprintf('\\phi_{Bottom} = %g', Phi_Bottom)}, ...
        'EdgeColor','none', 'FontSize', 12, 'FontWeight', 'bold');
    
    annotation('textbox', [0.91 0.52 0.08 0.05], ...
        'String', {sprintf('\\phi_{Right} = %g', Phi_Right)}, ...
        'EdgeColor','none', 'FontSize', 12, 'FontWeight', 'bold');

    %% Compute cell-center coordinates
    Px = 1:Nx;
    Py = 1:Ny;
    xP = (Px - 0.5) * dx;
    yP = (Py - 0.5) * dy;
    [Xp, Yp] = meshgrid(xP, yP);

    %% Extract internal phi values (remove ghost cells)
    phi_internal = phi(2:end-1, 2:end-1);

    %% 3D Surface Plot
    figure;
    surf(Xp, Yp, flipud(phi_internal), 'EdgeColor', 'none');
    colorbar;
    colormap(jet);
    title('3D Surface Plot of \phi at Cell Centers');
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('\phi (Temperature)');
    view(45, 30); 

    %% Display phi matrix
    fprintf('The Final T Matrix (internal and near-boundary nodes):\n');
    disp(phi_internal)
    disp('&');
    disp(phi)
    
    fprintf('Converged in iterations with error = %.5e\n', error);
end

