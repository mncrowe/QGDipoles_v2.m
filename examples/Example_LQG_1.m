% Plots a one-layer Larichev-Reznik dipole (LRD). This example appeared as
% Example 1 in:
%
% - Crowe & Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic 
%       model, J. Fluid Mech., 994, R1

addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;      % vortex speed
a = 1;      % vortex radius
R = 1;      % barotropic Rossby radius
beta = 1;   % background PV gradient

% Linear system parameters:

M = 8;     % set the number of coefficients calculated (default: 7)

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonLQG(Nz = 1, U = U, a = a, R_0 = R, beta = beta, M = M, ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, EVP_method = 'quad_EVP');

% Display K and plot psi and q:

disp(['K = ' num2str(modon.K)])

PlotModon(modon, 'psi', xlim = [-2 2], ylim = [-2 2], contour = 10)
PlotModon(modon, 'q', xlim = [-2 2], ylim = [-2 2], contour = 10)