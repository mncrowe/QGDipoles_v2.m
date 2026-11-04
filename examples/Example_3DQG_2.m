% Plots a barotropic modon in a 3D QG model and compares with the 
% equivalent one-layer solution. This example appeared as Example 2 in:
%
% - Crowe, 2026, Three-dimensional modons in a quasi-geostrophic model, J.
%       Geophys. Astrophys. Fluid Dyn.

addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;      % vortex speed
a = 1;      % vortex radius
H = 1;      % layer depth

% 3D modon parameters:

Nz = 21;

% Set K_0:

K_0_SQG = 4;
K_0_3DQG = [zeros(Nz-1, 1); 6];

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solutions:

modon_SQG = CreateModonSQG(U = U, a = a, H = H, Nx = Nx, Ny = Ny, ...
    Lx = Lx, Ly = Ly, K_0 = K_0_SQG);

modon_3DQG = CreateModon3DQG(Nz = Nz, U = U, a = a, Nx = Nx, Ny = Ny, ...
    Lx = Lx, Ly = Ly, K_0 = K_0_3DQG);

% Display K for the SQG model and K at z = 0 for the 3D QG model:

disp(['K_SQG = ' num2str(modon_SQG.K)])
disp(['K_3DQG(z = 0) = ' num2str(modon_3DQG.K(end))])

% Plot psi and b in the x-y plane for the SQG solution:

PlotModon(modon_SQG, 'psi', xlim = [-2 2], ylim = [-2 2], contour = 6)
PlotModon(modon_SQG, 'b', xlim = [-2 2], ylim = [-2 2], contour = 6)

% Plot psi and b in the x-y plane (z = 0) and y-z plane (x = 0) for the 3DQG solution:

PlotModon(modon_3DQG, 'psi', xlim = [-2 2], ylim = [-2 2], index = Nz, contour = 6)
PlotModon(modon_3DQG, 'b', xlim = [-2 2], ylim = [-2 2], index = Nz, contour = 6)

PlotModon(modon_3DQG, 'psi', ylim = [-2 2], plane = 'yz', contour = 6)
PlotModon(modon_3DQG, 'b', ylim = [-2 2], plane = 'yz', contour = 6)