% Matlab code for double-TM reconstruction
% Corresponding results are presented in FIG. 1 of the paper "Non-invasive 
% focusing and imaging in scattering media with a fluorescence-based
% transmission matrix"
% Authors: Antoine Boniface / Jonathan Dong / Sylvain Gigan 
% contact: antoine.boniface@lkb.ens.fr   //  version 02/2020
clear all, close all, clc


%% Parameters initialization

v.n_grains = 10;  % number of speckle grains
v.input_dim = 256;  % number of input pixels
v.grain_size = 3;  % size of the speckle grain
v.n_target = 4;  % number of targets  
v.img_size = v.grain_size * v.n_grains;  % size of the speckle image


%% Generation of speckle for every input mode

% Speckle Scat.1 (medium in transmission) & Scat.2 (same medium but in reflexion / epi detection)
T1 = generate_tm(v, v.input_dim);  % Scat. 1: T1 = TM from the SLM to CAM2 (control camera, in transmission)
T2 = generate_tm(v, v.n_target);  % Scat. 2: T2 = TM from the targets to CAM1 (epi-detection)

% Target positions
x = randi(round(0.95*size(T1,1)), v.n_target, 1)' + 1;
y = randi(round(0.8*size(T1,2)), v.n_target, 1)' + 1;

% Speckle fluo 1P emitted by incoherent targets, for the first mode of T1
I_cam_init = generate_fluo_img(v, T1(:, :, 1), T2, y, x);

% Display
% Excitation intensity
figure(1); subplot(121); 
imagesc(abs(T1(:, :, 1)).^2); 
title('|E_{exc}|^2(p)'); daspect([1 1 1]); axis off
hold on; plot(x(1:end),y(1:end),'w+','linewidth',3,'markersize',10);
% Fluorescence image
figure(1); subplot(122); 
imagesc(I_cam_init); 
title('I_{out}(p)'); daspect([1 1 1]); axis off


%% Random patterns on the SLM

% Generation of input/output data
display('Generation of input/output data')
v.n_patterns = 5 * v.input_dim;  % number of random patterns we send                                                              % size of the camera
I_out = zeros(v.n_patterns, v.img_size^2); 
% We send v.n_patterns random patterns on the SLM
for k=1:v.n_patterns
    new_phase = 2 * pi * rand(1,v.input_dim);
    slm_mask(k, :) = new_phase;  % we keep the SLM patterns for the phase retrieval
    excitation_field = zeros(v.img_size, v.img_size);
    % This part can be vectorized
    for j = 1:v.input_dim
        % coherent illumination speckle of the targets
        excitation_field = excitation_field + T1(:,:,j)*exp(1i * new_phase(j));
    end
    % generation of the low contrasted speckle pattern from the targets through the second scatt. medium
    I_cam = generate_fluo_img(v, excitation_field, T2, y, x);
    I_out(k, :) = reshape(I_cam, 1, v.img_size^2);
end


%% Non-negative Matrix Factorization (NMF)
display('Non-negative Matrix Factorization')
[W, H] = nnmf(I_out, v.n_target);


%% Phase Retrieval (PR)
display('Phase Retrieval')
n_iter = 5e2;  % default 5e2, usually high for simple gradient descent
step_size = 1e-6;  % default 1e-6, change if necessary
TM = zeros(v.n_target,v.input_dim);
for k=1:v.n_target
    y = W(:,k);
    y = y/norm(y);
    A = exp(1i * slm_mask);
    x_pr = phase_retrieval(y, A, n_iter, step_size);
    TM(k,:) = x_pr;
end


%% Phase Conjugation and Light Focusing on all the targets
display('Phase Conjugation')
sbr = zeros(v.n_target, 1);  % signal-to-background ratio
for jj=1:v.n_target
    E_target = zeros(v.n_target,1);
    E_target(jj) = 1;
    E_in = TM' * E_target;
    % Convert to phase only
    slm_phase = E_in ./ abs(E_in);
    excitation_field = zeros(v.img_size, v.img_size);
    for j = 1:v.input_dim
        excitation_field = excitation_field + T1(:,:,j) * slm_phase(j);
    end
    sbr(jj) = max(max(abs(excitation_field).^2))/mean2(abs(T1(:,:,1)).^2);
    % Display
    figure(13), subplot(1, v.n_target, jj),
    imagesc(abs(excitation_field).^2), axis off
    daspect([1 1 1]); title(['focus on target n° = ' num2str(jj)]);
end
