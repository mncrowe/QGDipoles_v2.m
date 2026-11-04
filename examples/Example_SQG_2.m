% Plots a steady, bottom intensified retrograde modon travelling along a
% topographic slope as described in:
%
% - Crowe & Johnson, 2023, The evolution of surface quasi-geostrophic 
%       modons on sloping topography, J. Fluid. Mech. 970, A10
%
% and plotted in Figures 1 and 2. Note that this paper used a different 
% form for K that the functions here. The corresponding value is K + mu, 
% where mu = lambda*a/U is the nondimensional topographic slope gradient.

addpath ..\functions\
addpath ..\

% Problem parameters:

U = 1;          % vortex speed
a = 1;          % vortex radius
lambda = 2;     % topographic slope

% Numerical plotting parameters:

Nx = 512;
Ny = 512;
Lx = 10;
Ly = 10;

% Get modon solution:

modon = CreateModonSQG(U = U, a = a, slope = lambda, location = 'bottom', ...
    Nx = Nx, Ny = Ny, Lx = Lx, Ly = Ly);

% Display K and plot b and psi:

disp(['K = ' num2str(modon.K - lambda * a / U)])

PlotModon(modon, 'b', xlim = [-2 2], ylim = [-2 2], contour = 6)
PlotModon(modon, 'psi', xlim = [-2 2], ylim = [-2 2], contours = 6)