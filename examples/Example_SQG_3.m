% Plots an SQG surface modon in a finite depth layer with a free surface
% and background PV gradient (beta-plane). We plot psi and Psi = psi + Uy.

addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;      % vortex speed
a = 1;      % vortex radius
g = 10;     % use finite gravity to set the barotropic Rossby number
beta = 1;   % background PV gradient
H = 1;      % layer depth

% Linear system parameters:

M = 10;     % set the number of coefficients calculated (default: 7)

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonSQG(U = U, a = a, g = g, beta = beta, H = H, ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, M = M);

% Display K and plot psi and psi + Uy:

disp(['K = ' num2str(modon.K)])

PlotModon(modon, 'psi', xlim = [-2 2], ylim = [-2 2], contours = 6)
PlotModon(modon, 'Psi', xlim = [-2 2], ylim = [-2 2], contours = 6)