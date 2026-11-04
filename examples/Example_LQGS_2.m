% Plots a 5-layer modon moving in the direction of a linear background 
% shear flow. This example appears as Example 2 in:
%
% - Crowe, 2026, N-Layer quasi-geostrophic modons in a vertically sheared 
%       flow, ArXiv pre-print

% note, choice of R and linear shear means only PV gradients in top and
% bottom layers

addpath ..\functions\
addpath ..\

% Problem parameters:

Nz = 4;                         % number of layers
U = 1;                          % vortex speed
a = 1;                          % vortex radius
R = [1 2 2 1];                  % baroclinic Rossby radii
beta = 0;                       % background PV gradient
U_n = [0.9 0.6 0.3 0];          % background flow

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonLQG(Nz = Nz, U = U, a = a, R = R, beta = beta, ...
    U_n = U_n, Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, save_linsys = true);

% Display K:

disp(['K = ' mat2str(modon.K)])

% Plot Psi(x, y) in each top and bottom layers and psi and q at x = 0:

Psi = modon.psi + modon.grid.y .* reshape((U - U_n), 1, 1, Nz);
psi_q_yi = cat(3, squeeze(modon.psi(Nx/2 + 1, :, :)), squeeze(modon.q(Nx/2 + 1, :, :)));

splot(Psi(:, :, [1 Nz]), modon.grid.x, modon.grid.y, circle_radius = 1, ...
    interpreter = 'latex', show_contours = [-0.9 -0.5 -0.2 -0.05 0.05 0.2 0.5 0.9], ...
    label_subplots = true, panel_array = [2 1], title = {'$\Psi_1(x, y)$', '$\Psi_4(x, y)$'}, ...
    xlim = [-1.5 1.5], ylim = [-1.5 1.5], clim = [-1 1], xlabel = '$x$', ylabel = '$y$', ...
    fig_size = [800 300])

splot(psi_q_yi, modon.grid.y, interpreter = 'latex', show_contours = 6, label_subplots = true, ...
    panel_array = [2 1], title = {'$\psi_i''(x = 0, y)$', '$q_i(x = 0, y)$'}, ...
    xlim = [-2 2], xlabel = '$y$', ylabel = '$i$', fig_size = [800 300], ...
    y_ticks = 1:Nz, flip_axis = [0 1], label_start = 3)