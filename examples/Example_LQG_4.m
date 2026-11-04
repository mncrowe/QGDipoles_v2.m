% Plots a N-layer modon with a range of parameters and active layers. This 
% example appeared as supplementary material with:
%
% - Crowe & Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic 
%       model, J. Fluid Mech., 994, R1

addpath ..\functions\
addpath ..\

% Problem parameters:

Nz = 5;                         % number of layers
U = 1;                          % vortex speed
a = 1;                          % vortex radius
R = [1 1 2 1 2];                % baroclinic Rossby radii
beta = [1 0 1 2 0];             % background PV gradient
active_layer = [1 1 0 1 1];     % isolated vortex core on middle layer

% Linear system parameters:

K_0 = [5 5 4 5 7].^2;   % initial guess for K
M = 8;                  % set the number of coefficients calculated (default: 7)

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 20;
Ly = 20;

% Get modon solution:

modon = CreateModonLQG(Nz = Nz, U = U, a = a, R = R, beta = beta, M = M, ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, active_layers = active_layer, ...
    K_0 = K_0);

% Display K and plot psi and q in each layer on individual plots:

disp(['K = ' mat2str(modon.K)])

PlotModon(modon, 'psi', xlim = [-2 2], ylim = [-2 2], contour = 10, individual = true)
PlotModon(modon, 'q', xlim = [-2 2], ylim = [-2 2], contour = 10, individual = true)