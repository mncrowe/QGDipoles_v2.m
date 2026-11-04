function [modon, grid, params] = CreateModon3DQG(params)
% Creates a modon solution to the 3D QG equations using the method 
% described in:
%
% - Crowe, 2026, Three-dimensional modons in a quasi-geostrophic model, J.
%       Geophys. Astrophys. Fluid Dyn.
%
% The 3D QG equations are:
%
% [∂_t - U ∂_x]q + J[ψ, q] + β ∂_x ψ = 0     for z ∈ (-H, 0),           [1]
%
% where ψ represents the vortex streamfunction and
%
% q = ∇²ψ + ∂_z[f²/N²(z) ∂_z ψ]              for z ∈ (-H, 0),
%
% is the potential vorticity anomaly. The streamfunction satifies the
% boundary conditions:
%
% [∂_t - U ∂_x](∂_z ψ + N²(0)/g ψ) + J[ψ, ∂_z ψ + N²(0)/g ψ] = 0 on z = 0,
% [∂_t - U ∂_x](∂_z ψ) + J[ψ, ∂_z ψ + N²(-H) S/f y] = 0          on z = -H.
%
% The top boundary condition allows for small height changes in the free
% surface and the bottom boundary condition incorporates a weak slope in
% the y direction. Throughout J[*,*] is the Jacobian derivative and the 
% problem parameters are:
%
%   - U: vortex speed (x direction)
%   - a: vortex radius
%   - β: background PV gradient
%   - f: Coriolis frequency
%   - N(z): vertical buoyancy frequency
%   - g: gravitational acceleration
%   - H: fluid depth
%   - S: bottom slope gradient
%
% Assuming steady vortex solutions, Equation [1] and the boundary 
% conditions can be written as:
%
%               J[ψ + Uy, q + βy] = 0     for z ∈ (-H, 0),
%    J[ψ + Uy, ∂_z ψ + N²(0)/g ψ] = 0      on z = 0,
% J[ψ + Uy, ∂_z ψ + N²(-H) S/f y] = 0      on z = -H,
%
% hence:
%
%               q + βy = F₀(ψ + Uy, z)   for z ∈ (-H, 0),             [2.a]
%    ∂_z ψ + N²(0)/g ψ = F₁(ψ + Uy)       on z = 0,                   [2.b]
% ∂_z ψ + N²(-H) S/f y = F₂(ψ + Uy)       on z = -H,                  [2.c]
%
% where:
%
% F₀(φ, z) = {β/U φ             r > a,
%            {-K²(z)/a² φ       r < a,
%
% F₁(φ) = {0                    r > a,
%         {K₁/H φ               r < a,
%
% F₂(φ) = {N²(-H)S/(fU) φ       r > a,
%         {K₂/H φ               r < a.
% 
% Equations [2.a]-[2.c] can be solved using a Hankel transform and an 
% expansion in terms of Zernike polynomials. Vertical derivatives are
% discretised through replacement with differentiation matrices
% corresponding to Chebyshev spectral methods or finite differences. The 
% value of [K₁ K²(z) K₂] is determined during the calculation. The solution 
% can be rotated through an angle α relative to the +ve x-axis and 
% translated to origin x₀. If Rossby waves exist with phase speed equal to 
% the vortex speed, steady vortex solutions do not exist.
%
% Function keyword arguments:
%
% - U: vortex speed in x direction (default: 1)
% - a: vortex radius in x-y plane (default: 1)
% - x_0: vortex center (default: [0 0])
% - alpha: vortex angle from x-axis (default: 0)
% - M: number of coefficients used to calculate solution (default: 7)
% - Nz: number of vertical gridpoints (default: 5)
% - beta: background PV gradient (default: 0)
% - N: vertical density gradient, N(z) (default: @(z) 1 + 0 * z)
% - H: layer depth (default: 1)
% - f: Coriolis parameter (default: 0)
% - g: gravitational acceleration (default: Inf)
% - slope: slope gradient S (default: 0)
% - spectral: use Chebyshev methods, not finite differences (default: true)
% - K_0: initial guess for K (default: [0 16 ... 16 0])
% - a_0: initial guess for coefficients (default: zeros(Nz, M))
% - Nx: number of x gridpoints (default: 256)
% - Ny: number of y gridpoints (default: 256)
% - Lx: length of x grid (default: 10)
% - Ly: length of y grid (default: 10)
% - int_tol: tolerance for numerical integration (default: 1e-6)
% - optim_tol: tolerance for solving eigenvalue problem (default: 1e-8)
% - AB_tol: tolerance for setting matrix entries to zero (default: 1e-8)
% - save_linsys: save linear system terms in output (default: false)
% - show_progress: display calculation progress (default: false)
%
% Function outputs:
%
% modon: structure containing modon solution
% grid: structure containing grid, included in modon if nargout < 3
% params: structure containing parameters, included in modon if nargout < 3
%
% -------------------------------------------------------------------------
% Note: The input parameter K_0 and output K refer to the estimated and
%       calculated vales of [K₁ K²(z_interior) K₂] where z_interior is the
%       vertical grid excluding the endpoints. Typically, K(z) ≈ 4 so 
%       taking K_0 ≈ [0 16 ... 16 0] is generally sufficient.
%
% Note: this function outputs ψ(x, y, z) and q(x, y, z) as defined above on 
%       the specified numerical grid. It also outputs b(x, y, z) = f ∂_z ψ.
%
% Note: for an LQG modon, index 3 of the output fields describes the layer
%       number, with 1 corresponding to the top layer and Nz to the bottom
%       layer. Conversely, for a 3D QG modon, index 3 of the output
%       corresponds to the value at height grid.z. Since grid.z is ordered
%       in increasing z, index 1 corresponds to the field value on z = -H,
%       i.e. the bottom, and index Nz corresponds to the field value on 
%       z = 0, i.e. the top.
% -------------------------------------------------------------------------

arguments
    params.U             (1,1) double             = 1
    params.a             (1,1) double             = 1
    params.x_0           (1,2) double             = [0, 0]
    params.alpha         (1,1) double             = 0
    params.M             (1,1) double {isinteger} = 7
    params.Nz            (1,1) double {isinteger} = 5
    params.beta          (1,1) double             = 0
    params.N             (1,1) function_handle    = @(z) 1 + 0*z;
    params.H             (1,1) double             = 1
    params.f             (1,1) double             = 1
    params.g             (1,1) double             = Inf
    params.slope         (1,1) double             = 0
    params.spectral      (1,1) logical            = true
    params.K_0           (:,1) double             = 16
    params.a_0           (:,:) double             = []
    params.Nx            (1,1) double {isinteger} = 256
    params.Ny            (1,1) double {isinteger} = 256
    params.Lx            (1,1) double             = 10
    params.Ly            (1,1) double             = 10
    params.int_tol       (1,1) double             = 1e-6
    params.optim_tol     (1,1) double             = 1e-8
    params.AB_tol        (1,1) double             = 1e-8
    params.save_linsys   (1,1) logical            = false
    params.show_progress (1,1) logical            = false
end

% set K_0 if not defined:

if numel(params.K_0) ~= params.Nz
    params.K_0 = params.K_0 *  [0 ones(1, params.Nz - 2) 0];
end

% Set a_0 = 0 if not specified:

if numel(params.a_0) == 0
    params.a_0 = zeros(params.Nz, params.M);
end

% Suppress singular matrix warning:

warn_state = warning;
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');

% Start timer:

tic;

% Extract parameters:

U = params.U;
a = params.a;
beta = params.beta;
alpha = params.alpha;

M = params.M;
Nz = params.Nz;

Nx = params.Nx;
Ny = params.Ny;

% Update progress:

if params.show_progress
    disp(' - Building linear system ...')
end

% Define parameters and functions:

lambda = @(zeta) a * params.f ./ (params.N(params.H * zeta) * params.H);
lambda_s = params.N(0) * sqrt(params.H / params.g);
gamma = params.slope * params.H * params.N(-params.H)^2 / (params.f * U);

mu = [-gamma beta*a^2/U*ones(1, Nz-2) 0];

if params.spectral

    [Dzeta, zeta] = grid_spectral(type = 'chebyshev2', N = Nz, L = [-1, 0]);

    Lambda = diag(lambda(zeta).^2);

    M0 = -Dzeta * Lambda * Dzeta;

else

    [Dzeta, zeta] = grid_finite_diff(N = Nz, L = [-1, 0], order = 1, stencil = [-1 0 1]);
    Dzeta_left = grid_finite_diff(N = Nz, L = [-1, 0], order = 1, stencil = [-1 0]);
    Dzeta_right = grid_finite_diff(N = Nz, L = [-1, 0], order = 1, stencil = [0 1]);

    Lambda = diag(lambda(zeta).^2);            

    M0 = -Dzeta_left * Lambda * Dzeta_right;

end

I0 = diag([0 ones(1, Nz-2) 0]);

delta = max(max(abs(M0))) / max(max(abs(Dzeta)));

M0(1, :) = delta * Dzeta(1, :);
M0(Nz, :) = delta * (Dzeta(Nz, :) + [zeros(1, Nz-1) lambda_s^2]);

mu(1) = delta * mu(1);

A_func = @(xi, I0, M0, C) (xi^2 * I0 + M0) / (xi^2 * I0 + M0 + C) / xi;
B_func = @(xi, I0, M0, C) inv(xi^2 * I0 + M0 + C) / xi;

% Build A and B matrices for linear system:

C = diag(mu);

A = zeros(Nz * M);
B0 = zeros(Nz * M);
iz = 1:Nz;

for j = 0:M-1

    for k = 0:j

        A_kj = JJ_int(@(x) A_func(x, I0, M0, C), j, k, tol = params.int_tol, atol = 1e-4 * params.int_tol);
        B_kj = JJ_int(@(x) B_func(x, I0, M0, C), j, k, tol = params.int_tol, atol = 1e-4 * params.int_tol);

        A(j*Nz + iz, k*Nz + iz) = A_kj;
        B0(j*Nz + iz, k*Nz + iz) = B_kj;

        A(k*Nz + iz, j*Nz + iz) = A_kj;
        B0(k*Nz + iz, j*Nz + iz) = B_kj;

    end

end

% Set very small matrix entries to zero:

A(abs(A) < params.AB_tol) = 0;
B0(abs(B0) < params.AB_tol) = 0;

% Construct remaining linear system terms:

B = cell(Nz, 1);
c = zeros(Nz * M, Nz + 1);
d = zeros(Nz * M, Nz);

c(1:Nz, 1) = 1/4 * mu;

for n = 1:Nz

    B{n} = kron(eye(M), diag(1:Nz == n)) * B0;
    c(:, n+1) = 1/4 * [(1:Nz) == n, zeros(1, (M-1)*Nz)];
    d(:, n) = kron((-1) .^ (0:M-1)', (1:Nz == n)');

end

% Update progress:

if params.show_progress
    disp(' - Solving linear system ...')
end

% Solve linear system for coefficients:

K_0 = params.K_0;
a_0 = reshape(params.a_0, Nz*M, 1);

K_0(1) = delta * K_0(1);
K_0(Nz) = delta * K_0(Nz);

[K, a_n, err] = EVP_optim(A, B, c, d, K_0 = K_0, a_0 = a_0, tol = params.optim_tol);

a_n = reshape(a_n, Nz, M);

K(1) = K(1) / delta;
K(Nz) = K(Nz) / delta;

% Save parameters and linear system solution in output structure:

modon.type = '3DQG';

if nargout < 3; modon.params = params; end

modon.K = K';
modon.a_n = a_n;
modon.err = rms(err);

% Update progress:

if params.show_progress
    disp(' - Calculating fields ...')
end

% Create numerical grid:

grid = create_grid([Nx Ny], [params.Lx params.Ly]);
grid.z = zeta * params.H;

r = sqrt((grid.x - params.x_0(1)).^2 + (grid.y - params.x_0(2)).^2);
theta = atan2(grid.y, grid.x);

% Calculate intermediate function, F, using coefficients, a_n:

F = zeros(Nx, Ny, Nz);

for j = 1:M
    for n = 1:Nz
        F(:, :, n) = F(:, :, n) + a_n(n, j) * R_n(r / a, j - 1);
    end
end

F = U / a * F .* sin(theta - alpha);

% Move to Fourier space:

Fh = FT_2D(F, 'forward');
Fh(Nx/2+1, Ny/2+1, :) = 0;     % set (0, 0) Fourier mode to zero

% Invert for psi then find b and q:

Delta = reshape(grid.k.^2 + grid.l.^2 + eps, 1, 1, Nx, Ny) .* I0 + reshape((M0 + C) / a^2, Nz, Nz, 1, 1);
Delta_inv = pageinv(Delta);
Delta_inv(isinf(Delta_inv)) = 0;

psih = pagemtimes(Delta_inv, reshape(permute(Fh, [3 1 2]), Nz, 1, Nx, Ny));
psi = real(FT_2D(permute(squeeze(psih), [2 3 1]), 'inverse'));

b = params.f * permute(squeeze(pagemtimes(Dzeta, reshape(permute(psi, [3 1 2]), Nz, 1, Nx, Ny))), [2 3 1]);

if params.spectral
    Delta_q = reshape(grid.k.^2 + grid.l.^2 + eps, 1, 1, Nx, Ny) .* eye(Nz) - ...
        reshape((Dzeta * Lambda * Dzeta) / a^2, Nz, Nz, 1, 1);
else
    Delta_q = reshape(grid.k.^2 + grid.l.^2 + eps, 1, 1, Nx, Ny) .* eye(Nz) - ...
        reshape((Dzeta_left * Lambda * Dzeta_right) / a^2, Nz, Nz, 1, 1);
end

qh = -pagemtimes(Delta_q, psih);
q = real(FT_2D(permute(squeeze(qh), [2 3 1]), 'inverse'));

% Save psi, b, q in output structure:

modon.psi = psi;
modon.b = b;
modon.q = q;

% Save linear system:

if params.save_linsys

    linsys.lambda = lambda;
    linsys.lambda_s = lambda_s;
    linsys.gamma = gamma;
    linsys.mu = mu;

    linsys.A = A;
    linsys.B = B;
    linsys.c = c;
    linsys.d = d;

    modon.linsys = linsys;
    
end

% Add runtime and grid to modon structure:

modon.runtime = toc;
if nargout == 1; modon.grid = grid; end

% Reactivate warnings:

warning(warn_state)

end