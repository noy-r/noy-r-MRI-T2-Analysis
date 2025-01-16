% Add the 'mapVBVD-main' folder to MATLAB's path %

addpath('mapVBVD-main');

% call mapVBVD directly 
twix = mapVBVD('meas_MID00172_FID189318_se.dat');

twix{:};

% In this tutorial we do not care for noise data, so let's get rid of the first measurement:
twix = twix{1};
% Now let's look at the object that holds the image information
twix.image;
%As you can see there is a lot of information, it seems that we acquired
%28 slices, one average and one segment. Here's the matrix size that would be created if we would call twix.image():
twix.image.dataSize(1:11);

%Let's get rid of the average dimension by telling the twix object to automatically average them during the
%read procedure. In addition, we also want to get rid of the oversampling that is almost always performed
%in read direction (2 x oversampling = 2 x larger FoV in read):

twix.image.flagDoAverage = true; % automatically average the average dimension
twix.image.flagRemoveOS = true;

%Now let's see how the two flags changed the 'size' of our data:
twix.image.dataSize(1:11)

%Now, let's look at different segments of our data. What are segments? Segments indicate how the acquisition
%of the data was 'segmented', e.g. in an EPI or TSE acquisition. This data is from a TSE and the
%number of segments is equal to the turbo factor. So each of this segments corresponds to a different echo
%in our TSE echo train.

data = twix.image( :, :, :, 1, 2, 1, 1, 1, 1, 1, :);
%Ah, young scribe, the discrepancy between your column number (256) and the data size column dimension (128) stems from the removal of oversampling in the readout direction.


%No Segments: There is no fourth dimension for data, confirming there’s only one segment.
%Only One Coil or Slice to Visualize: You’ll focus on visualizing data across the first three dimensions.

%% Single Coil Visualization - color
% single slice = 2
data = twix.image( :, :, :, 1, 2, 1, 1, 1, 1, 1, :);
figure();
coil = 1; % Select the coil to visualize
single_k_space = abs(squeeze(data(:, coil, :))); % Extract k-space magnitude for the coil

% Visualize k-space in color
% Exaggerate contrast with a gamma adjustment
% The signal is very low, thus powering by 1/x increasing it
imagesc(single_k_space.^0.2);
colorbar;
title(['Colored K-Space for Coil ', num2str(coil)]);
ylabel('Phase Encoding (Lin)');
xlabel('Readout (Col)');

%Visualize All Coils
%To visualize all 58 coils, create a montage:
figure();
num_coils = size(data, 2); % Total number of coils

sum_k_spaces=zeros(size(single_k_space));
for coil = 1:num_coils
    subplot(ceil(sqrt(num_coils)), ceil(sqrt(num_coils)), coil);
    k_space=abs(squeeze(data(:, coil, :))).^0.2;
    imagesc(k_space); % Show k-space for each coil 
    title(['Coil ', num2str(coil)]);
    axis off;
    sum_k_spaces=sum_k_spaces+k_space;
end
mean_k_spaces=sum_k_spaces./num_coils;


figure();
imagesc(mean_k_spaces);
colorbar;
title(['Mean K-spaces', num2str(coil)]);
ylabel('Phase Encoding (Lin)');
xlabel('Readout (Col)');

% The final matrix size of our image (single precision is enough):
os = twix.image.NCol/twix.image.dataSize(1); % oversampling factor
img = zeros([twix.image.NCol/os, twix.image.NLin, twix.image.NSli], 'single');
for sli = 1:twix.image.NSli
    % read in the data for slice 'sli'
    slice_data = twix.image(:,:,:,1,sli);
    % fft in col and lin dimension:
    fft_dims = [1 3];
    for f = fft_dims
        %Dimensions in the Data:
        %Dimension 1 (Col): Corresponds to the frequency encoding direction 
        %Dimension 3 (Lin): Corresponds to the phase encoding direction​
        slice_data = ifftshift(ifft(fftshift(slice_data,f),[],f),f);
        %fftshift: Centers the k-space data so that the zero-frequency component is in the middle.
        %ifft: Performs the Inverse FFT along the specified dimension.
        %ifftshift: Re-centers the data after the FFT operation.
    end
    % sum-of-square coil combination:
    % Coil data is combined into a single image for each slice using the sum-of-squares (SOS) method:
    img(:,:,sli) = squeeze(sqrt(sum(abs(slice_data).^2,2)));
end

% And here's the result for the three slices:
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4])
imagesc(flipud(img(:,:)), [0, 0.7*max(img(:))]), colormap gray,
% Purpose: Display the reconstructed image (img) as a scaled grayscale image.
% flipud(img(:,:)):
% Flips the image vertically (upside down to match typical medical imaging conventions).
% [0, 0.7*max(img(:))]:
% Scales the intensity range of the image.
% 0: The minimum intensity (black).
% 0.7*max(img(:)): Caps the maximum intensity to 70% of the image's brightest pixel, improving contrast by avoiding oversaturation.
% imagesc:Automatically scales the image to fit the axes.
axis image off;

TR = twix.hdr.MeasYaps.alTR{1}; %us
TEs = twix.hdr.MeasYaps.alTE; %us

sli_thickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
sli_position = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
turbo_factor = twix.hdr.MeasYaps.sFastImaging.lTurboFactor;

%% Calculate T2 map

% Initialize reconstructed images for TE1 and TE2
img_TE1 = zeros(size(img, 1), size(img, 2), 14, 'single'); % Slices 1–14
img_TE2 = zeros(size(img, 1), size(img, 2), 14, 'single'); % Slices 15–28

for sli = 1:twix.image.NSli
    % Read the data for each slice
    slice_data = twix.image(:,:,:,1,sli);
    
    % FFT along the frequency and phase encoding dimensions
    fft_dims = [1 3];
    for f = fft_dims
        slice_data = ifftshift(ifft(fftshift(slice_data, f), [], f), f);
    end
    
    % Combine coils using sum-of-squares (SOS)
    reconstructed_slice = squeeze(sqrt(sum(abs(slice_data).^2, 2)));
    
    % Assign to the appropriate TE group
    if sli <= 14
        img_TE1(:, :, sli) = reconstructed_slice;
    else
        img_TE2(:, :, sli - 14) = reconstructed_slice;
    end
end

% Convert TEs from microseconds to seconds
TE1 = TEs{1} / 1e6;
TE2 = TEs{2} / 1e6;

% Calculate T2 map
T2_map = (TE2 - TE1) ./ log(img_TE1 ./ img_TE2);

% Handle invalid values (e.g., division by zero)
T2_map(~isfinite(T2_map)) = 0;

% Visualize the T2 map
figure();
imagesc(mean(T2_map, 3)); % Average T2 map across slices
colormap('jet');
colorbar;
title('T2 Relaxation Time Map');
xlabel('X (mm)');
ylabel('Y (mm)');

% Visualize a single slice from TE1 and TE2
% slice_num = 7; % Select a slice to visualize
% 
% figure();
% subplot(1, 2, 1);
% imagesc(img_TE1(:, :, slice_num));
% colormap('gray');
% title(['Reconstructed Image (TE1, Slice ', num2str(slice_num), ')']);
% axis image off;
% 
% subplot(1, 2, 2);
% imagesc(img_TE2(:, :, slice_num));
% colormap('gray');
% title(['Reconstructed Image (TE2, Slice ', num2str(slice_num), ')']);
% axis image off;

% % Visualize T2 map for a single slice
% figure();
% imagesc(T2_map(:, :, slice_num));
% colormap('jet');
% colorbar;
% title(['T2 Map (Slice ', num2str(slice_num), ')']);
% xlabel('X (mm)');
% ylabel('Y (mm)');

%% k-space games
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
rec_img_TE1 = recostruct_from_k_space(kspace_TE1); % Reconstruct from original k-space (TE1)
rec_img_TE2 = recostruct_from_k_space(kspace_TE2); % Reconstruct from original k-space (TE2)

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
rec_img_TE1_LPF = recostruct_from_k_space(kspace_TE1_LPF); % Reconstruct from low-pass filtered k-space (TE1)
rec_img_TE2_LPF = recostruct_from_k_space(kspace_TE2_LPF); % Reconstruct from low-pass filtered k-space (TE2)

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
rec_img_TE1_HPF = recostruct_from_k_space(kspace_TE1_HPF); % Reconstruct from high-pass filtered k-space (TE1)
rec_img_TE2_HPF = recostruct_from_k_space(kspace_TE2_HPF); % Reconstruct from high-pass filtered k-space (TE2)

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

%% Visualizing as 3D Volume:
% % Define voxel dimensions% Define voxel dimensions
% voxel_size_x = 1; % Assuming isotropic voxels in x and y
% voxel_size_y = 1;
% voxel_size_z = sli_thickness; % Slice thickness in z (from the header)
% 
% % Scale axes
% [x, y, z] = meshgrid(1:size(img, 2), 1:size(img, 1), (0:size(img, 3)-1) * voxel_size_z);
% 
% % Visualize slices in the 3D volume
% figure();
% slice(x, y, z, img, ...
%     round(size(img, 2)/2), round(size(img, 1)/2), round(size(img, 3)/2)); % Slice at the middle
% colormap(gray);
% shading interp;
% title('3D Brain Volume (Interactive Slices)');
% xlabel('X (mm)');
% ylabel('Y (mm)');
% zlabel('Z (mm)');
% view(3); % 3D perspective
% axis tight;

function rec_img = recostruct_from_k_space(k_data) 
    fft_dims = [1 3];
        for f = fft_dims
            k_data = ifftshift(ifft(fftshift(k_data, f), [], f), f);
        end
    
        % Combine coils using sum-of-squares (SOS)
        rec_img = squeeze(sqrt(sum(abs(k_data).^2, 2)));
end
