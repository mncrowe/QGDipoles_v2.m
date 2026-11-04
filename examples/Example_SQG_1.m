% Plots an SQG modon solution at the surface of an infinite depth layer.
% These solutions appeared in:
%
% - Johnson & Crowe, 2023, Oceanic dipoles in a surface quasi-geostrophic 
%       model, J. Fluid. Mech. 958, R2

addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;  % vortex speed
a = 1;  % vortex radius

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonSQG(U = U, a = a, Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly);

% Display K and plot b and psi:

disp(['K = ' num2str(modon.K)])

PlotModon(modon, 'b', xlim = [-2 2], ylim = [-2 2])
PlotModon(modon, 'psi', xlim = [-2 2], ylim = [-2 2])