function [TM, cumul_TM] = generate_tm(v , n_speckle)
%GENSPECKLE generates a coherent Transmission Matrix
% It models the propagation of light through a complex scattering medium
% To obtain a finite speckle grain size (due to the diffraction limit)
% we use a low-pass filter in the Fourier space, associated with a pupil function

% Pupil definition
% This part may be vectorized
pupil = ones(v.n_grains, v.n_grains);
for j = 1:v.n_grains
    for k = 1:v.n_grains
        if (j-0.5-v.n_grains/2)^2 + (k-0.5-v.n_grains/2)^2 >= (v.n_grains/2)^2
            pupil(j, k) = 0;
        end
    end
end

% A speckle is modeled by a random phase in the Fourier space
bruit = exp(2*1i*pi * rand(v.n_grains, v.n_grains, n_speckle));
for j = 1:n_speckle
    bruit(:, :, j) = bruit(:, :, j) .* pupil;
end

% Fourier transform to go in the object space
% with zero padding for smoothing
TM = zeros(v.img_size, v.img_size, n_speckle);
for j = 1:n_speckle
    TM(:, :, j) = fft2(fftshift(padarray(bruit(:, :, j), ...
        [v.img_size/2 - v.n_grains/2, v.img_size/2 - v.n_grains/2])));
end

% I don't understand why we sum here TO CHECK
cumul_TM = zeros(v.img_size, v.img_size);
for j = 1:n_speckle
    cumul_TM = cumul_TM + TM(:, :, j);
end

end
