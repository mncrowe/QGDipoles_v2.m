function [K, a, err] = EVP_quad(A, B, c, d, options)
% Solves the problem (A-K*B)*a = c_0 + K*c_1 such that d^T*a = 0.
%
% A: M x M matrix
% B: M x M matrix
% c: 2 vectors (M x 1), enter as array c = [c_0 c_1]
% d: vectors (M x 1), enter as array d = [d]
%

arguments
    A             (:,:) double
    B             (:,:) double
    c             (:,2) double
    d             (:,1) double
    options.K_0   (1,1) double  = zeros(1, length(B))
end

% Build quadratic eigenvalue problem by inverting the better conditioned matrix:

if rcond(B) > rcond(A)
    M1 = (d' * (B \ c(:,1))) * A - (c(:,1) * d') * (B \ A);
    M2 = -(d' * (B \ c(:,1))) * B + (d' * (B \ c(:,2))) * A - (c(:,2) * d') * (B \ A);
    M3 = -(d' * (B \ c(:,2))) * B;
else
    M1 = (d' * (A \ c(:,1))) * A;
    M2 = (d' * (A \ c(:,2))) * A - (d' * (A \ c(:,1))) * B + (c(:,1) * d') * (A \ B);
    M3 = -(d' * (A \ c(:,2))) * B + (c(:,2) * d') * (A \ B);
end

% Put the system in canonical form and solve:

K = polyeig(M1, M2, M3);
K = K(~isinf(K));

% Find the eigenvalue closest to K_0:

[~, i] = min((K - options.K_0) .^ 2);

% Set K and a_n using the identified eigenvalue:

K = K(i);
a = (A - K * B) \ (c(:,1) + K * c(:,2));

% Estimate error as difference between left and right hand side:

err = rms((A - K * B) * a - c(:,1) - K * c(:,2));

end