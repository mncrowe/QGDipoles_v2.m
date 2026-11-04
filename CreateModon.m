function [modon, grid, params] = CreateModon(params)
% Creates a modon solution of the specified type.
%
% See function descriptions of CreateModonLQG.m, CreateModon3DQG.m and
% CreateModonSQG.m for description of each problem setup, method, and
% parameter list.
%
% Function keyword arguments:
%
% - type: 'LQG', '3DQG' or 'SQG' (default: 'LQG')
% - parameters: see corresponding function description for relevant
%               list of parameters.
%
% -------------------------------------------------------------------------
% Note: this function requires CreateModonLQG.m, CreateModon3DQG.m and
%       CreateModonSQG.m to create modons of the corresponding type.
% -------------------------------------------------------------------------

arguments
    params.type          (1,:) char               = 'LQG'
    params.U             (1,1) double             = 1
    params.a             (1,1) double             = 1
    params.x_0           (1,2) double             = [0, 0]
    params.alpha         (1,1) double             = 0
    params.M             (1,1) double {isinteger} = 7
    params.Nz            (1,1) double {isinteger} = 5
    params.R             (:,:) double             = 1
    params.beta          (1,:) double             = 0
    params.R_0           (1,1) double             = Inf
    params.U_n           (:,1) double             = 0
    params.calc_R        (1,1) logical            = false
    params.rho           (1,:) double             = []
    params.rho_0         (1,1) double             = 1025
    params.N             (1,1)                    = NaN
    params.H             (1,:) double             = 1
    params.f             (1,1) double             = 1
    params.g             (1,1) double             = Inf
    params.active_layers (:,1) logical            = true
    params.slope         (1,1) double             = 0
    params.location      (1,:) char               = 'surface'
    params.spectral      (1,1) logical            = true
    params.K_0           (:,1) double             = []
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

% Set N based on problem type if not specified:

if isnan(params.N)
    switch params.type
        case {'LQG', 'SQG'}
            params.N = 1;
        case '3DQG'
            params.N = @(z) 1 + 0*z;
    end
end

% Set default K_0 for SQG case:

if strcmp(params.type, 'SQG') && isempty(params.K_0)
    params.K_0 = 7;
end

% Call modon creation function based on type:

switch params.type

    case 'LQG'

        [modon, grid, params] = CreateModonLQG(U = params.U, a = params.a, x_0 = params.x_0, ...
            alpha = params.alpha, M = params.M, Nz = params.Nz, R = params.R, beta = params.beta, ...
            R_0 = params.R_0, U_n = params.U_n, K_0 = params.K_0, a_0 = params.a_0, Nx = params.Nx, ...
            Ny = params.Ny, Lx = params.Lx, Ly = params.Ly, int_tol = params.int_tol, ...
            optim_tol = params.optim_tol, AB_tol = params.AB_tol, save_linsys = params.save_linsys, ...
            calc_R = params.calc_R, f = params.f, g = params.g, H = params.H, rho = params.rho, ...
            rho_0 = params.rho_0, show_progress = params.show_progress, ...
            active_layers = params.active_layers, EVP_method = params.EVP_method);

    case 'SQG'

        [modon, grid, params] = CreateModonSQG(U = params.U, a = params.a, x_0 = params.x_0, ...
            alpha = params.alpha, M = params.M, beta = params.beta, N = params.N, H = params.H, ...
            f = params.f, g = params.g, K_0 = params.K_0, a_0 = params.a_0, Nx = params.Nx, ...
            Ny = params.Ny, Lx = params.Lx, Ly = params.Ly, int_tol = params.int_tol, ...
            optim_tol = params.optim_tol, AB_tol = params.AB_tol, save_linsys = params.save_linsys, ...
            show_progress = params.show_progress, slope = params.slope, location = params.location, ...
            EVP_method = params.EVP_method);

    case '3DQG'

        [modon, grid, params] = CreateModon3DQG(U = params.U, a = params.a, x_0 = params.x_0, ...
            alpha = params.alpha, M = params.M, Nz = params.Nz, beta = params.beta, N = params.N, ...
            H = params.H, f = params.f, g = params.g, slope = params.slope, spectral = params.spectral, ...
            K_0 = params.K_0, a_0 = params.a_0, Nx = params.Nx, Ny = params.Ny, Lx = params.Lx, ...
            Ly = params.Ly, int_tol = params.int_tol, optim_tol = params.optim_tol, ...
            AB_tol = params.AB_tol, save_linsys = params.save_linsys, ...
            show_progress = params.show_progress);
        
end

% Add grid and params to output structure if not included as output arguments:

if nargout < 2
    modon.grid = grid;
end

if nargout < 3
    modon.params = params;
end

end