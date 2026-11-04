function [K, a, err] = EVP_optim(A, B, c, d, options)
% Solves the problem (A-sum[K_i*B_i])*a = c_0 + sum[K_i*c_i] such that
% d_k^T*a = 0 for i,k = {1,2,...,N}.
%
% A: M x M matrix
% B: cell array of N M x M matrices (B = cell(N,1); B{i} = B_i)
% c: vectors (M x 1), enter as array c = [c_0 c_1 c_2 .. c_N]
% d: vectors (M x 1), enter as array d = [d_1 d_2 .. d_N]
%
% keyword options:
% K_0: initial guess for [K_1, K_2, ..., K_N]
% a_0: initial guess for a, vector of length M
% debug: shows iteration output if set to true (default: false)
% tol: function tolerance for fsolve (default: 1e-10)
%
% This function works by defining:
% 
%       f(x) = (A-sum[K_i*B_i])*a - c_0 - sum[K_i*c_i],
%
% where x = [K_1, K_2, .. K_N, a'_1, .. a'_{M-N}] and the a'_i are
% coefficients of the M-N basis vectors which span the space orthogonal to
% [d_1, d_2, ..., d_N]. We now find the solution to f(x) = 0  using the
% Matlab fsolve function from the Optimisation tooolbox. The solution x
% gives the values of K_i and allows a to be reconstructed from the a'_i. 
% Note that the complexity of this problem does not increase as the number
% of eigenvalues, N, does since a new eigenvalue introduces a new condition
% which reduces the dimension of the space containing a. Resonable guesses
% for K_0 are generally sufficient, allowing a_0 to be left as 0.

arguments
    A             (:,:) double
    B             (:,1) cell
    c             (:,:) double
    d             (:,:) double
    options.K_0   (1,:) double  = zeros(1, length(B))
    options.a_0   (:,1) double  = zeros(size(A, 1), 1)
    options.debug (1,1) logical = false
    options.tol   (1,1) double  = 1e-10
end

% first find a basis for the space orthogonal to the d_j:

[e, V, I] = orthog_space(d);
[N, J] = size(e);

% create initial guess randomly or from supplied data if given:

x0 = V \ options.a_0;
x0 = [reshape(options.K_0, N-J, 1); x0(I)];

% find solution x which solves f(x) = 0:

f = @(x) EVP_optim_f(x, A, B, c, e);

if options.debug
    optim_options = optimoptions('fsolve', TolFun = options.tol / N, SpecifyObjectiveGradient = true, TolX = options.tol / N, display = 'iter');
else
    optim_options = optimoptions('fsolve', TolFun = options.tol / N, SpecifyObjectiveGradient = true, TolX = options.tol / N, display = 'none');
end

[x, err] = fsolve(f, x0, optim_options);

% determine K and a from x:

K = reshape(x(1:N-J), 1, N-J);
a = e * x(N-J+1:N);

end

% definition of the function f(x):

function [f_out, Jac] = EVP_optim_f(x, A, B, c, e)
       
    [N, J] = size(e);

    a = e*x(N-J+1:N);
    M = A;
    v = c(:, 1);
    Jac = zeros(N, N-J);

    for i = 1:N-J
        M = M - x(i) * B{i};
        v = v + x(i) * c(:, i+1);
        Jac(:, i) = -B{i} * a - c(:, i+1);
    end

    f_out = (M * a - v);
    Jac = [Jac M * e];
    
end