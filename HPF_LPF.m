clc , clear, close all

%% Add the 'mapVBVD-main' folder to MATLAB's path %
addpath('mapVBVD-main');

% call mapVBVD directly 
twix = mapVBVD('meas_MID00172_FID189318_se.dat');
% In this tutorial we do not care for noise data, so let's get rid of the first measurement:
twix = twix{1};
twix.image.flagDoAverage = true; % automatically average the average dimension 
twix.image.flagRemoveOS = true; %% automatically removes oversampling

%% K space games 

slice_num = 7; % Select a slice to visualize
kspace_TE1 = twix.image(:, :, :, 1, slice_num);  % Extract for TE1 (slice 7)
kspace_TE2 = twix.image(:, :, :, 1, slice_num + slice_num*2); % Extract for TE2 (slice 21)

% Define k-space grid (for filtering purposes)
[rows,coils, cols] = size(kspace_TE1(:, :, :));
[X, Y] = meshgrid(-cols/2:(cols/2 - 1), -rows/2:(rows/2 - 1));

% Define filter radii (you can modify these for your needs)
radius_LPF = rows / 6; % For Low-Pass Filter (1/4th of k-space)
radius_HPF = rows / 6; % For High-Pass Filter (smaller radius for edges)

% Create LPF and HPF masks
LPF_mask = sqrt(X.^2 + Y.^2) <= radius_LPF; % Circular low-pass filter mask
HPF_mask = sqrt(X.^2 + Y.^2) > radius_HPF;  % Circular high-pass filter mask

% Ensure the filter mask is 128x1x128
LPF_mask_resized = reshape(LPF_mask, [128, 1, 128]); % Reshape to 128x1x128
HPF_mask_resized = reshape(HPF_mask, [128, 1, 128]); % Reshape to 128x1x128


% % Resize the filters to match k-space size (128x128) — coils will stay the same
% LPF_mask_resized = imresize(LPF_mask, [rows, cols]);
% HPF_mask_resized = imresize(HPF_mask, [rows, cols]);

figure();
subplot(1, 2, 1);
imagesc(squeeze(LPF_mask_resized)); 
colormap('gray'); 
title('Low-Pass Filter (LPF) Mask');
axis image off;

subplot(1, 2, 2);
imagesc(squeeze(HPF_mask_resized));
colormap('gray');
title('High-Pass Filter (HPF) Mask');
axis image off;

% This figure shows the LPF and HPF masks as separate images.

kspace_TE1_LPF = kspace_TE1;  % Start with original k-space
kspace_TE1_HPF = kspace_TE1;  % Start with original k-space

kspace_TE2_LPF = kspace_TE2;  % Start with original k-space
kspace_TE2_HPF = kspace_TE2;  % Start with original k-space

for coil = 1:coils
    % Apply the filter to each coil separately
    kspace_TE1_LPF(:, coil, :) = kspace_TE1(:, coil, :) .* LPF_mask_resized;  % Apply LPF to each coil of TE1
    kspace_TE1_HPF(:, coil, :) = kspace_TE1(:, coil, :) .* HPF_mask_resized;  % Apply HPF to each coil of TE1

    kspace_TE2_LPF(:, coil, :) = kspace_TE2(:, coil, :) .* LPF_mask_resized;  % Apply LPF to each coil of TE2
    kspace_TE2_HPF(:, coil, :) = kspace_TE2(:, coil, :) .* HPF_mask_resized;  % Apply HPF to each coil of TE2
end
% Visualize k-space with and without filters for TE1 and TE2
figure();
subplot(3, 2, 1);
imagesc(squeeze(abs(kspace_TE1(:, 1, :))).^0.2); % Original k-space (TE1) - showing only the first coil for example
title('Original k-space (TE1)');
axis image off;
colormap('gray');

subplot(3, 2, 2);
imagesc(squeeze(abs(kspace_TE2(:, 1, :))).^0.2); % Original k-space (TE2) - showing only the first coil for example
title('Original k-space (TE2)');
axis image off;
colormap('gray');

subplot(3, 2, 3);
imagesc(squeeze(abs(kspace_TE1_LPF(:, 1, :))).^0.2); % Low-pass filtered k-space (TE1)
title('Low-Pass Filtered k-space (TE1)');
axis image off;
colormap('gray');

subplot(3, 2, 4);
imagesc(squeeze(abs(kspace_TE2_LPF(:, 1, :))).^0.2); % Low-pass filtered k-space (TE1)
title('Low-Pass Filtered k-space (TE2)');
axis image off;
colormap('gray');

subplot(3, 2, 5);
imagesc(squeeze(abs(kspace_TE1_HPF(:, 1, :))).^0.2); % High-pass filtered k-space (TE1)
title('High-Pass Filtered k-space (TE1)');
axis image off;
colormap('gray');

subplot(3, 2, 6);
imagesc(squeeze(abs(kspace_TE2_HPF(:, 1, :))).^0.2); % Low-pass filtered k-space (TE1)
title('High-Pass Filtered k-space (TE2)');
axis image off;
colormap('gray');

% Reconstruct the image from original k-space (TE1 and TE2)
rec_img_TE1 = RecostructFromKSpaceSiemens(kspace_TE1); % Reconstruct from original k-space (TE1)
rec_img_TE2 = RecostructFromKSpaceSiemens(kspace_TE2); % Reconstruct from original k-space (TE2)

% Visualize reconstructed images from original k-space
figure();
subplot(3, 2, 1);
imagesc(rec_img_TE1);
title('TE1, Original');
axis image off;
colormap('gray');

subplot(3, 2, 2);
imagesc(rec_img_TE2);
title('TE2, Original');
axis image off;
colormap('gray');

% Reconstruct the image from low-pass filtered k-space (TE1 and TE2)
rec_img_TE1_LPF = RecostructFromKSpaceSiemens(kspace_TE1_LPF); % Reconstruct from low-pass filtered k-space (TE1)
rec_img_TE2_LPF = RecostructFromKSpaceSiemens(kspace_TE2_LPF); % Reconstruct from low-pass filtered k-space (TE2)

% Visualize reconstructed images from low-pass filtered k-space
subplot(3, 2, 3);
imagesc(rec_img_TE1_LPF);
title('TE1, LPF');
axis image off;
colormap('gray');

subplot(3, 2, 4);
imagesc(rec_img_TE2_LPF);
title('TE2, LPF');
axis image off;
colormap('gray');

% Reconstruct the image from high-pass filtered k-space (TE1 and TE2)
rec_img_TE1_HPF = RecostructFromKSpaceSiemens(kspace_TE1_HPF); % Reconstruct from high-pass filtered k-space (TE1)
rec_img_TE2_HPF = RecostructFromKSpaceSiemens(kspace_TE2_HPF); % Reconstruct from high-pass filtered k-space (TE2)

% Visualize reconstructed images from high-pass filtered k-space
subplot(3, 2, 5);
imagesc(rec_img_TE1_HPF);
title('TE1, HPF');
axis image off;
colormap('gray');

subplot(3, 2, 6);
imagesc(rec_img_TE2_HPF);
title('TE2, HPF)');
axis image off;
colormap('gray');