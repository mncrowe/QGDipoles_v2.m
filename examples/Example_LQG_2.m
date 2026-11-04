% Plots a two-layer baroclinic modon on a beta-plane with a bottom slope. 
% This example appeared as Example 2 in:
%
% - Crowe & Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic 
%       model, J. Fluid Mech., 994, R1

addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;          % vortex speed
a = 1;          % vortex radius
R = [1 1];      % baroclinic Rossby radii
beta = [0 1];   % background PV gradient

% Linear system parameters:

M = 12;         % set the number of coefficients calculated (default: 7)

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonLQG(Nz = 2, U = U, a = a, R = R, beta = beta, M = M, ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly);

% Display K and plot psi and q in each layer in the same window:

disp(['K = (' num2str(modon.K(1)) ', ' num2str(modon.K(2)) ')'])

PlotModon(modon, 'psi', xlim = [-2 2], ylim = [-2 2], contour = 10)
PlotModon(modon, 'q', xlim = [-2 2], ylim = [-2 2], contour = 10)