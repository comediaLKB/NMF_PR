function I_cam = generate_fluo_img(v, Speckleinit, Speckle_collec, Y, X)
%GENSPECKLEFLUO generates the fluorescence images collected by the camera
% It uses the previously-defined Transmission Matrices in Speckle[...]
% It sums incoherently the fluorescent speckles of the different targets

% Excitation of each target
for bi=1:size(Y,2)
	Fluo(bi) = abs(Speckleinit(Y(bi),X(bi))).^2;
end

% Fluorescence image of each target
for p = 1:v.n_target
    I_collec(:,:,p) = abs(Speckle_collec(:,:,p)).^2;
    I_fluo(:,:,p) = I_collec(:,:,p)*Fluo(p);
end

% Summing the fluorescence
I_cam = sum(I_fluo,3);

end
