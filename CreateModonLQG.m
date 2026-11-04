function [modon, grid, params] = CreateModonLQG(params)
% Creates a modon solution to the N-layer QG equations using the method 
% described in:
%
% - Crowe & Johnson, 2024, Modon solutions in an N-layer quasi-geostrophic 
%       model, J. Fluid Mech., 994, R1
% - Crowe, 2025, QGDipoles.jl: A Julia package for calculating dipolar 
%       vortex solutions to the Quasi-Geostrophic equations, JOSS, 10 (108)
% - Crowe, 2026, N-Layer quasi-geostrophic modons in a vertically sheared 
%       flow, ArXiv pre-print
%
% The N-layer QG equations are:
%
% [∂_t - U ∂_x]qᵢ + J[ψᵢ - Uᵢy, qᵢ - ΛᵢⱼUⱼy] + βᵢ ∂_x ψᵢ = 0,           [1]
%
% where ψᵢ represents the vortex streamfunction and
%
% qᵢ = ∇²ψᵢ + Λᵢⱼψⱼ,
% 
% is the potential vorticity anomaly. Here J[*,*] is the Jacobian 
% derivative and the problem parameters are:
%
%   - U: vortex speed (x direction)
%   - a: vortex radius
%   - Uᵢ: background velocity (x direction) in layer i
%   - βᵢ: background PV gradient in layer i
%   - Λᵢⱼ: matrix of Rossby radii coupling layers i and j
% 
% The matrix Λ is tridiagonal and given by:
%
% Λ = [-[R₀⁻² + R₁.₅,₁⁻²]       R₁.₅,₁⁻²                0            ...  ]
%     [    R₁.₅,₂⁻²     -[R₁.₅,₂⁻² + R₂.₅,₂⁻²]       R₂.₅,₂⁻²        ...  ]
%     [       0                 R₂.₅,₃⁻²     -[R₂.₅,₃⁻² + R₃.₅,₃⁻²]  ...  ]
%     [      ...                  ...                  ...           ...  ]
%
% Where Rᵢ₋₀.₅,ⱼ = √[g(ρᵢ - ρᵢ₋₁)Hⱼ/(ρ₀f²)] is the barotropic radius
% describing the effect on the flow in layer j resulting from a deformation
% of the interface between layers i-1 and i. Similarly, R₀ = √[gH₁/f²] is
% the barotropic Rossby radius describing the effects of the deformation of
% the top free surface on the flow in the top layer (layer 1).
%
% Assuming steady vortex solutions, Equation [1] can be written as:
%
% J[ψᵢ + (U - Uᵢ)y, qᵢ + (βᵢ - ΛᵢⱼUⱼ)y] = 0,
%
% Hence:
%
% qᵢ + (βᵢ - ΛᵢⱼUⱼ)y = Fᵢ(ψᵢ + (U - Uᵢ)y),                              [2]
%
% where we take ψᵢ + (U - Uᵢ)y = 0 on the vortex boundary and
%
% Fᵢ(z) = {-Kᵢ²/a² z                :   r < a
%         {(βᵢ - ΛᵢⱼUⱼ)/(U - Uᵢ) z  :   r > a
%
% Equation [2] can be solved using a Hankel transform and an expansion in
% terms of Zernike polynomials. The value of Kᵢ² is determined during the
% calculation. Passive layers in which Kᵢ² = -(βᵢ - ΛᵢⱼUⱼ)a²/(U - Uᵢ) and
% there are no isolated regions of vorticity in layer(s) i can also be
% found. The solution can be rotated through an angle α relative to the +ve
% x-axis and translated to origin x₀. If Rossby waves exist with phase
% speed equal to the vortex speed, steady vortex solutions do not exist.
%
% Function keyword arguments:
%
% - U: vortex speed in x direction (default: 1)
% - a: vortex radius in x-y plane (default: 1)
% - x_0: vortex center (default: [0 0])
% - alpha: vortex angle from x-axis (default: 0)
% - M: number of coefficients used to calculate solution (default: 7)
% - Nz: number of vertical layers (default: 3)
% - R: baroclinic Rossby radii [Rᵢ₋₀.₅,ᵢ; Rᵢ₊₀.₅,ᵢ] (default: ones(2, Nz))
% - beta: background PV gradient in each layer (default: [0 ... 0])
% - R_0: barotropic Rossby radius in top layer (default: Inf)
% - U_n: backgrouynd velocity in x direction (default: [0 ... 0])
% - calc_R: calculate R and R_0 using physical parameters (default: false)
% - rho: density in each layer (default: UNSPECIFIED)  [calc_R = true only]
% - rho_0: average layer density (default: 1025)       [calc_R = true only]
% - H: layer depth (default: UNSPECIFIED)              [calc_R = true only]
% - f: Coriolis parameter (default: 1e-4)              [calc_R = true only]
% - g: gravitational acceleration (default: 9.8)       [calc_R = true only]
% - active_layers: layer has an active vortex (default: [true ... true])
% - K_0: initial guess for K (default: 20)
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
% - EVP_method: method for solving EVP, 'quad_EVP' or 'optim' (default)
%
% Function outputs:
%
% modon: structure containing modon solution
% grid: structure containing grid, included in modon if nargout < 3
% params: structure containing parameters, included in modon if nargout < 3
%
% -------------------------------------------------------------------------
% Note: The input parameter K_0 and output K refer to the estimated and
%       calculated vales of Kᵢ², not Kᵢ. Typically, Kᵢ ≈ 4-5 so taking 
%       K_0 ≈ 20 is generally sufficient.
%
% Note: this function outputs ψᵢ(x, y) and qᵢ(x, y) as defined above on the
%       specified numerical grid.
%
% Note: the underlying eigenvalue problem may be solved by two methods; an
%       optimisation approach and conversion to a quadratic eigenvalue
%       problem which may be solved directly. The methods are EVP_method = 
%       'quad_EVP' or EVP_method = 'optim'. The 'quad_EVP' method is 
%       generally faster and more accurate but can only be used for modons
%       with one active layer. The default is set as 'optim' for
%       generality. Different values of K_0 may be needed for each method.
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
    params.Nz            (1,1) double {isinteger} = 3
    params.R             (:,:) double             = 1
    params.beta          (1,:) double             = 0
    params.R_0           (1,1) double             = Inf
    params.U_n           (:,1) double             = 0
    params.calc_R        (1,1) logical            = false
    params.rho           (1,:) double             = []
    params.rho_0         (1,1) double             = 1025
    params.H             (1,:) double             = []
    params.f             (1,1) double             = 1e-4
    params.g             (1,1) double             = 9.8
    params.active_layers (:,1) logical            = true
    params.K_0           (:,1) double             = 20
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
    params.EVP_method    (1,:) char               = 'optim'
end

% Set a_0 = 0 if not specified:

if numel(params.a_0) == 0
    params.a_0 = zeros(params.Nz, params.M);
end

% Set R_0 and R using physical parameters if calc_R = true:

if params.calc_R
    if params.Nz > 1
        if numel(params.H) > 0 && numel(params.rho) > 0
            params.R_0 = sqrt(params.g * params.H(1)) ./ params.f;
            params.R = sqrt(params.g * params.H .* [params.rho - params.rho([1 1:end-1]); ...
                params.rho([2:end end]) - params.rho] / params.rho_0) / params.f;
        else
            error('Specify H and rho to calculate Rossby radii when calc_R = true.')
        end
    else
        if numel(params.H) == 1
            params.R_0 = sqrt(params.g * params.H(1)) ./ params.f;
        else
            error('Specify H to calculate Rossby radius when calc_R = true.')
        end
    end
end

% Create array for K_0, R, and vectors for beta and U_n if vectors/scalars are entered:

if numel(params.K_0) ~= params.Nz
    params.K_0 = params.K_0 * ones(1, params.Nz);
end

switch numel(params.R)
    case params.Nz
        params.R = [reshape(params.R, 1, []); reshape(params.R, 1, [])];
    case 2 * params.Nz
    otherwise
        params.R = params.R(1) * ones(2, params.Nz);
end

if length(params.beta) ~= params.Nz
    params.beta = params.beta(1) * ones(1, params.Nz);
end

if length(params.U_n) ~= params.Nz
    params.U_n = params.U_n(1) * ones(params.Nz, 1);
end

if length(params.active_layers) ~= params.Nz
    params.active_layers = repmat(params.active_layers, params.Nz, 1);
end

params.R(1, 1) = NaN;
params.R(2, end) = NaN;

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
U_n = params.U_n;

M = params.M;
Nz = params.Nz;

Nx = params.Nx;
Ny = params.Ny;

% Update progress:

if params.show_progress
    disp(' - Building linear system ...')
end

% Define parameters and functions:

R_ij = -diag([0 params.R(1, 2:end).^-2] + [params.R(2, 1:end-1).^-2 0]) + ...
        diag(params.R(2, 1:end-1).^-2, 1) + diag(params.R(1, 2:end).^-2, -1);

R_ij(1, 1) = R_ij(1, 1) - 1 / params.R_0^2;

M0 = -a^2 * R_ij;
mu = (beta' - R_ij * U_n) * a^2 ./ (U - U_n);
gamma = (U - U_n) / U;

I = eye(Nz);

A_func = @(xi, M0, C) (xi^2 * I + M0) / (xi^2 * I + M0 + C) / xi;
B_func = @(xi, M0, C) inv(xi^2 * I + M0 + C) / xi;

% Build A and B matrices for linear system:

C = diag(mu);

A = zeros(Nz * M);
B0 = zeros(Nz * M);
iz = 1:Nz;

for j = 0:M-1

    for k = 0:j

        A_kj = JJ_int(@(x) A_func(x, M0, C), j, k, tol = params.int_tol, atol = 1e-4 * params.int_tol);
        B_kj = JJ_int(@(x) B_func(x, M0, C), j, k, tol = params.int_tol, atol = 1e-4 * params.int_tol);

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

B = zeros(Nz * M, Nz * M, Nz);
c = zeros(Nz * M, Nz + 1);
d = zeros(Nz * M, Nz);

c(1:Nz, 1) = 1/4 * mu .* gamma;

for n = 1:Nz

    B(:, :, n) = kron(eye(M), diag(1:Nz == n)) * B0;
    c(:, n+1) = 1/4 * gamma(n) * [(1:Nz) == n, zeros(1, (M-1)*Nz)];
    d(:, n) = kron((-1) .^ (0:M-1)', (1:Nz == n)');

end

% Remove rows and columns corresponding to passive layers:

if any(~params.active_layers)

    i_1 = logical(kron(ones(M, 1), params.active_layers));
    
    A = A(i_1, i_1);
    B = B(i_1, i_1, params.active_layers);
    c = c(i_1, [true; params.active_layers]);
    d = d(i_1, params.active_layers);

    params.a_0 = params.a_0(params.active_layers, :);
    params.K_0 = params.K_0(params.active_layers);

end

% Convert B to cell array:

B = num2cell(B, [1 2]);

% Update progress:

if params.show_progress
    disp(' - Solving linear system ...')
end

% Solve linear system for coefficients:

a_0 = reshape(params.a_0, (Nz - sum(~params.active_layers)) * M, 1);

switch params.EVP_method
    case 'quad_EVP'
        if length(params.K_0) == 1
            [K_p, a_n_p, err] = EVP_quad(A, B{1}, c, d, K_0 = params.K_0);
        else
            error('EVP_method = ''quad_EVP'' requires only 1 active layer')
        end
    case 'optim'
        [K_p, a_n_p, err] = EVP_optim(A, B, c, d, K_0 = params.K_0, a_0 = a_0, tol = params.optim_tol);
    otherwise
        error('EVP_method not recognised.')
end

a_n_p = reshape(a_n_p, Nz - sum(~params.active_layers), M);

% Include passive layers in solution:

K = -mu;
K(params.active_layers) = K_p;

a_n = zeros(Nz, M);
a_n(params.active_layers, :) = a_n_p;

% Raise warning if number of passive layers in solution doesn't match input:

if sum(K + mu < 10 * params.optim_tol) ~= sum(~params.active_layers)
    warning('The number of active layers in the solution does not match the specified number of active layers.')
end

% Save parameters and linear system solution in output structure:

modon.type = 'LQG';

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

if Nz == 1; F = reshape(F, Nx, Ny, 1); end

% Move to Fourier space:

F_FT = FT_2D(F,'forward');
F_FT(Nx/2+1, Ny/2+1, :) = 0;     % set (0, 0) Fourier mode to zero

% Invert for psi then find q:

Delta = reshape(grid.k.^2 + grid.l.^2 + eps, 1, 1, Nx, Ny) .* I + reshape((M0 + C) / a^2, Nz, Nz, 1, 1);
Delta_inv = pageinv(Delta);
Delta_inv(isinf(Delta_inv)) = 0;

psih = pagemtimes(Delta_inv, reshape(permute(F_FT, [3 1 2]), Nz, 1, Nx, Ny));
psi = real(FT_2D(squeeze(permute(psih, [3 4 1 2])), 'inverse'));

Delta_q = reshape(grid.k.^2 + grid.l.^2 + eps, 1, 1, Nx, Ny) .* I + reshape(M0 / a^2, Nz, Nz, 1, 1);
    
qh = -pagemtimes(Delta_q, psih);

q = real(FT_2D(squeeze(permute(qh, [3 4 1 2])), 'inverse'));

% Save psi, q in output structure:

modon.psi = psi;
modon.q = q;

% Save linear system:

if params.save_linsys

    linsys.mu = mu;
    linsys.R_ij = R_ij;
    linsys.gamma = gamma;

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