% Plots a baroclinic beta-plane modon on a topographic slope with variable 
% stratification. This example appeared as Example 3 in:
%
% - Crowe, 2026, Three-dimensional modons in a quasi-geostrophic model, J.
%       Geophys. Astrophys. Fluid Dyn.

addpath ..\functions\
addpath ..\

% Problem parameters:

Nz = 21;                % number of z gridpoints
U = 1;                  % vortex speed
a = 1;                  % vortex radius
H = 1;                  % layer depth
N = @(z) exp(z/2);      % depth-dependent stratification
beta = 1;               % background PV gradient
slope = 1;              % slope gradient in y direction

% Set K_0 guess:

K_0 = [-3; -beta/U*ones(Nz-2, 1); 5];

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModon3DQG(Nz = Nz, U = U, a = a, H = H, N = N, beta = beta, ...
    slope = slope, Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly, K_0 = K_0);

% Display K on z = -1 and z = 0:

disp(['K(-1) = ' num2str(modon.K(1))])
disp(['K(0) = ' num2str(modon.K(Nz))])

% Plot psi and b in the x-y plane (z = -1) and y-z plane (x = 0) for the 3DQG solution:

PlotModon(modon, 'psi', xlim = [-3 3], ylim = [-3 3], index = 1, contour = [-2.2 -1.2 -0.2 0 0.2 1.2 2.2])
PlotModon(modon, 'b', xlim = [-3 3], ylim = [-3 3], index = 1, contour = [-2 -0.15 0.15 2])

PlotModon(modon, 'psi', ylim = [-3 3], plane = 'yz', contour = 7)
PlotModon(modon, 'b', ylim = [-3 3], plane = 'yz', contour = [-2 -0.15 0.15 2])