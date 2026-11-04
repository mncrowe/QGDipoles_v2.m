% Plots a barotropic modon in a 3D QG model and compares with the 
% equivalent one-layer solution. This example appeared as Example 1 in:
%
% - Crowe, 2026, Three-dimensional modons in a quasi-geostrophic model, J.
%       Geophys. Astrophys. Fluid Dyn.

%close all; clear
addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;
a = 1;
beta = 1;

% 3D modon parameters:

N = @(z) 10 + 0 * z;
Nz = 21;

% 1-layer modon parameters:

R = Inf;

% Set K_0 and a_0 guesses:

K_0_2D = 16;
K_0_3D = [0; 15.38*ones(Nz-2, 1); 0];
a_0_3D = [-12.55 -16.18 -4.026 -0.4213 -0.0244 -0.000905 0.0000335] .* [0; ones(Nz-2, 1); 0];

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 20;
Ly = 20;

% Get modon solutions:

modon_2D = CreateModonLQG(Nz = 1, U = U, a = a, R_0 = R, beta = beta, ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, K_0 = K_0_2D);

modon_3D = CreateModon3DQG(Nz = Nz, U = U, a = a, N = N, beta = beta, ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, K_0 = K_0_3D, a_0 = a_0_3D);

% Display K and plot 2D psi in the x-y plane:

disp(['K_{2D} = ' mat2str(modon_2D.K)])
disp(['mean(K_{3D}) = ' mat2str(mean(modon_3D.K(2:end-1)))])

PlotModon(modon_2D, 'psi', xlim = [-2 2], ylim = [-2 2], contour = 6)

% Plot 3D psi and q in the x-y plane at z = -1 and z = 0 and in the y-z plane along x = 0:

PlotModon(modon_3D, 'psi', xlim = [-2 2], ylim = [-2 2], index = [1 Nz], contour = 6)
PlotModon(modon_3D, 'q', xlim = [-2 2], ylim = [-2 2], index = [1 Nz], contour = 6)

PlotModon(modon_3D, 'psi', ylim = [-2 2], plane = 'yz', contour = 6)
PlotModon(modon_3D, 'q', ylim = [-2 2], plane = 'yz', contour = 6)