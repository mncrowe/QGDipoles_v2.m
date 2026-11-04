function grid = create_grid(N, L)
% Creates periodic Fourier grid centered on the origin
%
% N: [Nx Ny], number of gridpoints (default: [128 128])
% L: [Lx Ly], domain size (default: [10 10])

arguments
    N (1,2) double {isinteger} = [128 128]
    L (1,2) double             = [10 10]
end

Nx = N(1); Ny = N(2);
Lx = L(1); Ly = L(2);

dx = Lx / Nx;
dy = Ly / Ny;

grid.x = (-Lx/2:dx:(Lx/2-dx))';
grid.y = -Ly/2:dy:(Ly/2-dy);

grid.k = (pi/dx*(-1:2/Nx:(1-2/Nx)))';
grid.l = pi/dy*(-1:2/Ny:(1-2/Ny));

end