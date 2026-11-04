function [modon, grid, params] = CreateModonSQG(params)
% Creates a modon solution to the surface QG equations using the method 
% described in:
%
% - Johnson & Crowe, 2023, Oceanic dipoles in a surface quasi-geostrophic 
%       model, J. Fluid. Mech. 958, R2
% - Crowe & Johnson, 2023, The evolution of surface quasi-geostrophic 
%       modons on sloping topography, J. Fluid. Mech. 970, A10
% - Crowe, 2025, QGDipoles.jl: A Julia package for calculating dipolar
%       vortex solutions to the Quasi-Geostrophic equations, JOSS, 10 (108)
%
% The method presented here generalises and combines the methods presented
% in the work cited above.
%
% The 3D QG equations are:
%
% [∂_t - U ∂_x]q + J[ψ, q] + β ∂_x ψ = 0     for z ∈ (H₁, H₂),          [1]
%
% where ψ represents the vortex streamfunction and
%
% q = ∇²ψ + f²/N² ∂_z² ψ                     for z ∈ (H₁, H₂),
%
% is the potential vorticity anomaly. The streamfunction satifies the
% boundary conditions:
%
% [∂_t - U ∂_x](∂_z ψ + N²/g ψ) + J[ψ, ∂_z ψ + N²/g ψ] = 0     on z = H₁,
% [∂_t - U ∂_x](∂_z ψ) + J[ψ, ∂_z ψ + N² S/f y] = 0            on z = H₂.
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
%   - N: vertical buoyancy frequency
%   - g: gravitational acceleration
%   - H: fluid depth
%   - S: bottom slope gradient
%
% We now look for steady modon solutions so Equation [1] and the boundary
% conditions become:
%
%          J[ψ + Uy, q + βy] = 0     for z ∈ (H₁, H₂),
%  J[ψ + Uy, ∂_z ψ + N²/g ψ] = 0      on z = H₁,
% J[ψ + Uy, ∂_z ψ + N²S/f y] = 0      on z = H₂,
%
% hence:
%
%          q + βy = F₀(ψ + Uy)       for z ∈ (H₁, H₂),                [2.a]
%  ∂_z ψ + N²/g ψ = F₁(ψ + Uy)        on z = H₁,                      [2.b]
% ∂_z ψ + N²S/f y = F₂(ψ + Uy)        on z = H₂.                      [2.c]
%
% We now impose the SQG assumption:
%
% F₀(φ) = β/U φ,
%
% hence the interior PV arises only from the deflection around the vortex 
% in the y direction. Equation [2.a] reduces to the usual SQG assumption of 
% zero interior PV in the case of β = 0. We have
%
% ∇²ψ + f²/N² ∂_z² ψ - β/U ψ = 0       for z ∈ (H₁, H₂)                 [3]
%
% which can be solved subject to a boundary condition on the top or bottom.
% The other boundary will contain the vortex core.
%
% We now consider two case; surface modons, and bottom modons.
%
% Case 1: Surface modons:
%
% Here we assume an isolated region of PV on the top surface. We take:
%
% F₁(φ) = {NK₁/(fa) φ   r < a,
%         {0            r > a,
%
% F₂(φ) = N²S/(fU) φ,
%
% and note that for this case we have (H₁, H₂) = (0, -H) where H may be ∞.
% The appropriare boundary condition for Equation [3] comes from F₂ as:
%
% ∂_z ψ - N²S/(fU) ψ = 0      on z = -H.                                [4]
%
% We're left with a single equation on the top surface, z = 0, given by
% Equation [2.b]:
%
% ∂_z ψ + N²/g ψ = {NK₁/(fa) (ψ + Uy)  r < a,                           [5]
%                  {0                  r > a,
%
% where ∂_z ψ and ψ on z = 0 are related by a Dirichlet-to-Neumann operator
% which can be found by solving [3] using the boundary condition given
% in Equation [4].
%
% Case 2: Bottom modons:
%
% Here we assume an isolated region of PV on the bottom surface. We take:
%
% F₁(φ) = 0,
%
% F₂(φ) = {[N²S/(fU) - NK₂/(fa)] φ   r < a,
%         {N²S/(fU) φ                 r > a,
%
% and note that for this case we have (H₁, H₂) = (H, 0) where H may be ∞.
% The appropriate boundary condition for Equation [3] comes from F₁ as:
%
% ∂_z ψ + N²/g ψ = 0          on z = H.                                 [6]
%
% We're left with a single equation on the bottom surface, z = 0, given by
% Equation [2.c]:
%
% ∂_z ψ + N²S/f y = {[N²S/(fU) - NK₂/(fa)] (ψ + Uy)   r < a,            [7]
%                   {N²S/(fU) (ψ + Uy)                r > a,
%
% where ∂_z ψ and ψ on z = 0 are related by a Dirichlet-to-Neumann operator
% which can be found by solving [3] using the boundary condition given
% in Equation [6].
%
% For Cases 1 and 2, the Equations [5] and [7] are solved respectively. We
% use a Hankel transform and an expansion in terms of Zernike polynomials. 
% The value of K₁ or K₂ is determined during the calculation. The solution 
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
% - beta: background PV gradient (default: 0)
% - N: vertical density gradient, N(z) (default: 1)
% - H: layer depth (default: Inf)
% - f: Coriolis parameter (default: Inf)
% - g: gravitational acceleration (default: Inf)
% - slope: slope gradient S (default: 0)
% - location: Case from above; 'surface' or 'bottom' (default: 'surface')
% - K_0: initial guess for K (default: ±7 depending on 'location')
% - a_0: initial guess for coefficients (default: zeros(M, 1))
% - Nx: number of x gridpoints (default: 256)
% - Ny: number of y gridpoints (default: 256)
% - Lx: length of x grid (default: 10)
% - Ly: length of y grid (default: 10)
% - int_tol: tolerance for numerical integration (default: 1e-6)
% - optim_tol: tolerance for solving eigenvalue problem (default: 1e-8)
% - AB_tol: tolerance for setting matrix entries to zero (default: 1e-8)
% - save_linsys: save linear system terms in output (default: false)
% - show_progress: display calculation progress (default: false)
% - EVP_method: method for solving EVP, 'quad_EVP' (default) or 'optim'
%
% Function outputs:
%
% modon: structure containing modon solution
% grid: structure containing grid, included in modon if nargout < 3
% params: structure containing parameters, included in modon if nargout < 3
%
% -------------------------------------------------------------------------
% Note: The input parameter K_0 and output K refer to the estimated and
%       calculated vales of K₁ or K₂. Typically, K₁,K₂ ≈ 5 so taking 
%       K_0 ≈ 5 is generally sufficient.
%
% Note: this function outputs ψ(x, y, z = 0). It also outputs 
%       b(x, y, z = 0) = f ∂_z ψ and q(x, y, z = 0) = β/U ψ.
%
% Note: the underlying eigenvalue problem may be solved by two methods; an
%       optimisation approach and conversion to a quadratic eigenvalue
%       problem which may be solved directly. The methods are EVP_method = 
%       'quad_EVP' or EVP_method = 'optim'. For SQG solutions, the
%       'quad_EVP' method is generally faster and more accurate. Different
%       values of K_0 may be needed for each method.
% -------------------------------------------------------------------------

arguments
    params.U             (1,1) double             = 1
    params.a             (1,1) double             = 1
    params.x_0           (1,2) double             = [0, 0]
    params.alpha         (1,1) double             = 0
    params.M             (1,1) double {isinteger} = 7
    params.beta          (1,1) double             = 0
    params.N             (1,1) double             = 1
    params.H             (1,1) double             = Inf
    params.f             (1,1) double             = 1
    params.g             (1,1) double             = Inf
    params.slope         (1,1) double             = 0
    params.location      (1,:) char               = 'surface'
    params.K_0           (1,1) double             = 5
    params.a_0           (:,1) double             = []
    params.Nx            (1,1) double {isinteger} = 256
    params.Ny            (1,1) double {isinteger} = 256
    params.Lx            (1,1) double             = 10
    params.Ly            (1,1) double             = 10
    params.int_tol       (1,1) double             = 1e-6
    params.optim_tol     (1,1) double             = 1e-8
    params.AB_tol        (1,1) double             = 1e-8
    params.save_linsys   (1,1) logical            = false
    params.show_progress (1,1) logical            = false
    params.EVP_method    (1,:) char               = 'quad_EVP'
end

% Set a_0 = 0 if not specified:

if numel(params.a_0) == 0
    params.a_0 = zeros(params.M, 1);
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

Nx = params.Nx;
Ny = params.Ny;

% Update progress:

if params.show_progress
    disp(' - Building linear system ...')
end

% Define parameters and functions:

lambda = a * params.f / (params.N * params.H);       % a/R (inverse baroclinic Rossby radius)
lambda_p = a * params.N * params.f / params.g;       % a/R' (inverse barotropic Rossby radius)
mu = beta * a^2 / U;                                 % dimensionless background PV gradient
gamma = a * params.N * params.slope / U;             % dimensionless slope gradient

switch params.location
    case 'surface'
        C = [mu, lambda, gamma, lambda_p];
    case 'bottom'
        C = [mu, lambda, lambda_p, gamma];
    otherwise
        error('Modon location must be ''surface'' or ''bottom''.')
end

% Get Dirichlet-to-Neumann operator, simplify if possible to reduce computation time:

if C(2) == 0
    if C(1) == 0
        D = @(xi, C) xi;
    else
        D = @(xi, C) sqrt(xi.^2 + C(1));
    end
else
    if C(3) == 0
        if C(1) == 0
            D = @(xi, C) xi .* tanh((eps + xi) / C(2));
        else
            D = @(xi, C) sqrt(xi.^2 + C(1)) .* tanh(sqrt(eps + xi.^2 + C(1)) / C(2));
        end
    else
        if C(1) == 0
            D = @(xi, C) xi .* (C(3) + xi .* tanh((eps + xi) / C(2))) ./ (C(3) * tanh((eps + xi) / C(2)) + eps + xi);
        else
            D = @(xi, C) sqrt(xi.^2 + C(1)) .* (C(3) + sqrt(xi.^2 + C(1)) .* tanh(sqrt(eps + xi.^2 + C(1)) / C(2))) ./ ...
                    (C(3) * tanh(sqrt(eps + xi.^2 + C(1)) / C(2)) + sqrt(eps + xi.^2 + C(1)));
        end
    end
end

if ~all(C([1 2 4]) == 0)
    B_func = @(xi, C) 1 ./ ((D(xi, C(1:3)) + C(4)) .* xi);
end

% Build A and B matrices for linear system:

B0 = zeros(M);
A = diag(1 ./ (1:M) / 4);

for j = 0:M-1

    for k = 0:j

        if all(C([1 2 4]) == 0)
            B_kj = 4 * (-1)^(j-k+1) / ((2*j-2*k-1) * (2*j-2*k+1) * (2*j+2*k+3) * (2*j+2*k+5)) / pi;
        else
            B_kj = JJ_int(@(x) B_func(x, C), j, k, tol = params.int_tol, atol = 1e-4 * params.int_tol);
        end

        B0(j + 1, k + 1) = B_kj;
        B0(k + 1, j + 1) = B_kj;

    end

end

% Set very small matrix entries to zero:

A(abs(A) < params.AB_tol) = 0;
B0(abs(B0) < params.AB_tol) = 0;

% Construct remaining linear system terms:

B = cell(1);
B{1} = B0;

c = zeros(M, 2);
c(1, 2) = 1/4;

d = (-1).^(0:M-1)';

% Update progress:

if params.show_progress
    disp(' - Solving linear system ...')
end

% Solve linear system for coefficients:

switch params.EVP_method
    case 'quad_EVP'
        [K, a_n, err] = EVP_quad(A, B{1}, c, d, K_0 = params.K_0);
    case 'optim'
        [K, a_n, err] = EVP_optim(A, B, c, d, K_0 = params.K_0, a_0 = params.a_0, tol = params.optim_tol);
    otherwise
        error('EVP_method not recognised.')
end

a_n = reshape(a_n, 1, M);

% Save parameters and linear system solution in output structure:

modon.type = 'SQG';

params.Nz = 1;
if nargout < 3; modon.params = params; end

modon.K = K;
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

F = zeros(Nx, Ny);

for j = 1:M
    F(:, :) = F(:, :) + a_n(j) * R_n(r / a, j - 1);
end

F = a * U * F .* sin(theta - alpha);

% Move to Fourier space:

Fh = FT_2D(F, 'forward');
Fh(Nx/2+1, Ny/2+1, :) = 0;     % set (0, 0) Fourier mode to zero

% Invert for psi then find b and q:

dz = D(a * sqrt(grid.k.^2 + grid.l.^2), C(1:3));    % dimensionless vertical derivative

psih = Fh ./ (dz + C(4) + eps);
psi = real(FT_2D(psih, 'inverse'));

bh = params.N * dz .* psih / a;                     % factor of a fixes dimensions as dz dimensionless
b = FT_2D(bh, 'inverse');

if strcmp(params.location, 'bottom'); b = -b; end   % flip sign of b on bottom boundary as d/dz switches sign
        
q = beta / U * psi;

% Save psi, b, q in output structure:

modon.psi = psi;
modon.b = b;
modon.q = q;

% Save linear system:

if params.save_linsys

    linsys.mu = mu;
    linsys.lambda = lambda;
    linsys.lamdba_p = lambda_p;
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