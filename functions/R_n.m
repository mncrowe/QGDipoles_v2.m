function y = R_n(x, n)
% Calculates R_n(x), the Zernike radial function of order n, using Jacobi
% polynomials of order n

arguments
    x       double
    n (1,1) double {mustBeInteger}
end

y = (-1) ^ n * x .* jacobiP2(n, 0, 1, 2 * x.^2 - 1) .* (x <= 1);

end