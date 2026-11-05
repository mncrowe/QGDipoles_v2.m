% Plots a 3-layer modon in a depth-independent backgroynd flow with a 
% corresponding free-surface-induced PV gradient in the upper layer. This
% example appears as Example 1 in: 
%
% - Crowe, 2026, N-Layer quasi-geostrophic modons in a vertically sheared 
%       flow, ArXiv pre-print

% note constant shear still interesting for finite R_0 as there's a
% corresponding height change leading to a PV gradient

addpath ..\functions\
addpath ..\

% Problem parameters:

Nz = 3;                         % number of layers
U = 1;                          % vortex speed
a = 1;                          % vortex radius
R_0 = 1;                        % barotropic Rossby radius
R = [1 1 1];                    % baroclinic Rossby radii
beta = [1 1 0];                 % background PV gradient
U_n = [0.5 0.5 0.5];            % background flow

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonLQG(Nz = Nz, U = U, a = a, R_0 = R_0, R = R, beta = beta, ...
    U_n = U_n, Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, save_linsys = true);

% Display K:

disp(['K = ' mat2str(modon.K)])

% Calculate PV gradient:

disp(['β_i - Λ_{ij} U_j  = ' mat2str(modon.params.beta' - modon.linsys.R_ij * modon.params.U_n)])

% Plot psi and Psi in each layer:

F = zeros(Nx, Ny, 2*Nz);
F(:, :, 1:3) = modon.psi;
F(:, :, 4:6) = modon.psi + (U - U_n(1)) * modon.grid.y;

splot(F, modon.grid.x, modon.grid.y, circle_radius = 1, interior_labels = false, ...
    interpreter = 'latex', show_contours = [-0.75 -0.5 -0.25 0.25 0.5 0.75], label_subplots = true, ...
    panel_array = [3 2], title = {'$\psi_1''(x, y)$', '$\psi_2''(x, y)$', ...
    '$\psi_3''(x, y)$', '$\Psi_1(x, y)$', '$\Psi_2(x, y)$', '$\Psi_3(x, y)$'}, ...
    xlim = [-2 2], ylim = [-2 2], clim = [-1 1], xlabel = '$x$', ylabel = '$y$', ...
    fig_size = [1000 450])