% Plots a 6-layer modon moving against a quadratic background shear flow on
% a beta-plane. This example appears as Example 3 in:
%
% - Crowe, 2026, N-Layer quasi-geostrophic modons in a vertically sheared 
%       flow, ArXiv pre-print

% note choice of R and nonlinear shear gives PV gradients in middle layers,
% look at mu

addpath ..\functions\
addpath ..\

% Problem parameters:

Nz = 6;                             % number of layers
U = 1;                              % vortex speed
a = 1;                              % vortex radius
R = [1 2 3 3 2 1; 1 3 4 4 3 1];     % baroclinic Rossby radii
beta = 1;                           % background PV gradient
U_n = -[1 0.8 0.6 0.4 0.2 0].^2;    % background flow

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Set K_0:

K_0 = 30;

% Get modon solution:

modon = CreateModonLQG(Nz = Nz, U = U, a = a, R = R, beta = beta, ...
    U_n = U_n, Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, save_linsys = true, K_0 = K_0);

% Display K:

disp(['K = ' mat2str(modon.K)])

% Calculate PV gradient and R:

disp(['β_i - Λ_{ij} U_j  = ' mat2str(modon.params.beta' - modon.linsys.R_ij * modon.params.U_n)])
R_vec = reshape(modon.params.R, [], 1)';
disp(['R = ' mat2str(R_vec(2:end-1))])

% Plot Psi(x, y) in top and bottom layers and psi and q at x = 0:

Psi = modon.psi + modon.grid.y .* reshape((U - U_n), 1, 1, Nz);
psi_q_yi = cat(3, squeeze(modon.psi(Nx/2 + 1, :, :)), squeeze(modon.q(Nx/2 + 1, :, :)));

splot(Psi(:, :, [1 Nz]), modon.grid.x, modon.grid.y, circle_radius = 1, ...
    interpreter = 'latex', show_contours = [-2 -1.5 -1 -0.5 0.5 1 1.5 2], ...
    label_subplots = true, panel_array = [2 1], title = {'$\Psi_1(x, y)$', '$\Psi_6(x, y)$'}, ...
    xlim = [-1.5 1.5], ylim = [-1.5 1.5], clim = [-2 2], xlabel = '$x$', ylabel = '$y$', ...
    fig_size = [800 300])

splot(psi_q_yi, modon.grid.y, interpreter = 'latex', show_contours = 6, label_subplots = true, ...
    panel_array = [2 1], title = {'$\psi_i''(x = 0, y)$', '$q''_i(x = 0, y)$'}, ...
    xlim = [-2 2], xlabel = '$y$', ylabel = '$i$', fig_size = [800 300], ...
    y_ticks = 1:Nz, flip_axis = [0 1], label_start = 3)