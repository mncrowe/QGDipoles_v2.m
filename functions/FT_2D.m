function F = FT_2D(F, dir)
% Performs a 2D Fourier transform across the first two dimensions of F
%
% F: multidimensional array
% dir: direction; 'forward' or 'inverse'

arguments
    F   (:,:,:) double
    dir (1,:)   char   = 'forward'
end

switch dir
    case 'forward'
        F = fftshift(fftshift(fft2(F), 1), 2);
    case 'inverse'
        F = ifft2(fftshift(fftshift(F, 1), 2));
end

end