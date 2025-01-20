clc , clear, close all

% This code is combine all the tutorial maethods Noy and Rom did 
% other scripts are the division of the methods each to single script 

%% Add the 'mapVBVD-main' folder to MATLAB's path %
addpath('mapVBVD-main');

% call mapVBVD directly 
twix = mapVBVD('meas_MID00172_FID189318_se.dat');
% In this tutorial we do not care for noise data, so let's get rid of the first measurement:
twix = twix{1};
twix.image.flagDoAverage = true; % automatically average the average dimension 
twix.image.flagRemoveOS = true; %% automatically removes oversampling


%% K space Visualization - color
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
    rec_img = RecostructFromKSpaceSiemens(slice_data);
    img(:,:,sli)=rec_img;
end

% And here's the result for the three slices:
scrsz = get(groot,'ScreenSize');
figure('Position',[1 scrsz(4)/4 scrsz(3)/2 scrsz(4)/4])
imagesc(flipud(img(:,:)), [0, 0.7*max(img(:))]), colormap gray,
axis image off;




%% T2
TR = twix.hdr.MeasYaps.alTR{1}; %us
TEs = twix.hdr.MeasYaps.alTE; %us

sli_thickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;
sli_position = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
turbo_factor = twix.hdr.MeasYaps.sFastImaging.lTurboFactor;


% Pre-allocate arrays for TE1 and TE2 reconstructed images.
num_slices_TE = 14;
[cols, rows] = deal(twix.image.dataSize(1), twix.image.dataSize(3));  % Adjust based on your dataSize
img_TE1 = zeros(cols, rows, num_slices_TE, 'single');
img_TE2 = zeros(cols, rows, num_slices_TE, 'single');

for sli = 1:twix.image.NSli
    % Extract k-space data for the current slice
    slice_kspace = twix.image(:,:,:,1,sli);
    
    % Use your reconstruction function
    reconstructed_slice = RecostructFromKSpaceSiemens(slice_kspace);
    
    % Assign slice to TE1 or TE2 group based on slice index
    if sli <= num_slices_TE
        img_TE1(:,:,sli) = reconstructed_slice;
    elseif sli <= 2*num_slices_TE
        img_TE2(:,:,sli - num_slices_TE) = reconstructed_slice;
    else
        warning('Unexpected slice index beyond expected TE groups.');
    end
end

% Retrieve Echo Times and Compute T2 Map
% Retrieve Echo Times from header
if numel(TEs) < 2
    error('Not enough TE values found in header.');
end

TE1 = TEs{1} / 1e6;  % Convert to seconds
TE2 = TEs{2} / 1e6;

% Pre-allocate T2_map array for all slices (14 slices assumed for T2 map)
T2_map = zeros(cols, rows, num_slices_TE, 'single');

% Loop over each slice and compute pixel-wise T2
for slice_idx = 1:num_slices_TE
    S1 = img_TE1(:,:,slice_idx);
    S2 = img_TE2(:,:,slice_idx);
    
    % Apply a threshold to avoid division by zero issues
    min_signal = 1e-6;
    S1(S1 < min_signal) = min_signal;
    S2(S2 < min_signal) = min_signal;
    
    ratio = S1 ./ S2;
    valid = ratio > 0;  % Only compute where ratio > 0
    
    t2_slice = zeros(size(S1), 'single');
    t2_slice(valid) = (TE2 - TE1) ./ log(ratio(valid));
    t2_slice(~isfinite(t2_slice)) = 0;
    
    T2_map(:,:,slice_idx) = t2_slice;
end

% % Visualize T2 Map for each slice
%  Visualizing 14 figures
% for slice_idx = 1:num_slices_TE
%     figure;
%     imagesc(T2_map(:,:,slice_idx));
%     colormap('jet');
%     colorbar;
%     title(sprintf('T2 Map for Slice %d', slice_idx));
%     axis image off;
%     pause(0.5);
% end

% Determine subplot grid dimensions based on number of slices for T2 maps
num_plots = num_slices_TE; 
rows_subplot = ceil(sqrt(num_plots));
cols_subplot = ceil(num_plots / rows_subplot);

figure;
for slice_idx = 1:num_slices_TE
    subplot(rows_subplot, cols_subplot, slice_idx);
    imagesc(T2_map(:,:,slice_idx));
    colormap('jet');        % Use 'jet' colormap for T2 maps
    colorbar;               % Add a colorbar to each subplot
    title(sprintf('T2 Slice %d', slice_idx));
    axis image off;
end
sgtitle('T2 Maps for Slices 1–14');  % Adds a common title for the figure (available in newer MATLAB versions)

% Visualize mean T2 across slices if desired
mean_T2_across_slices = mean(T2_map, 3);
figure;
imagesc(mean_T2_across_slices);
colormap('jet');
colorbar;
title('Mean T2 Map Across Slices 1–14');
axis image off;

%% Relative Proton Density (PD) Mapping based on T2 Map and TE1 Images

% Pre-allocate PD_map array (relative PD) for all slices corresponding to TE1
PD_map = zeros(cols, rows, num_slices_TE, 'single');

for slice_idx = 1:num_slices_TE
    % Get signal at TE1 for current slice
    S1 = img_TE1(:,:,slice_idx);
    
    % Get corresponding T2 values for the slice
    T2_slice = T2_map(:,:,slice_idx);
    
    % Avoid division by zero or very small T2 values
    min_t2 = 1e-6;
    T2_slice(T2_slice < min_t2) = min_t2;
    
    % Compute relative PD using the equation:
    % PD_rel(x,y) ≈ S(x,y; TE1) * exp(TE1 / T2(x,y))
    PD_map(:,:,slice_idx) = S1 .* exp(TE1 ./ T2_slice);
end


% % Visualize PD Map for each slice
%  Visualizing 14 figures
% for slice_idx = 1:num_slices_TE
%     figure;
%     imagesc(PD_map(:,:,slice_idx));
%     colormap('gray');
%     colorbar;
%     title(sprintf('Relative PD Map for Slice %d', slice_idx));
%     axis image off;
%     pause(0.5);
% end

% Determine subplot grid dimensions based on number of slices
num_plots = num_slices_TE; 
rows_subplot = ceil(sqrt(num_plots));
cols_subplot = ceil(num_plots / rows_subplot);

figure;
for slice_idx = 1:num_slices_TE
    subplot(rows_subplot, cols_subplot, slice_idx);
    imagesc(PD_map(:,:,slice_idx));
    colormap('gray');
    colorbar;
    title(sprintf('Rel PD Slice %d', slice_idx));
    axis image off;
end
sgtitle('Relative PD Maps for Slices 1–14');  % Adds a common title for the figure

% Optionally, visualize mean PD across slices
mean_PD_across_slices = mean(PD_map, 3);
figure;
imagesc(mean_PD_across_slices);
colormap('gray');
colorbar;
title('Mean Relative PD Map Across Slices 1–14');
axis image off;


%% Define the slice index for TE1 you want to plot
slice_idx = 7;

% Create a new figure
figure;

% Display the image for slice 7 at TE1
imagesc(img_TE1(:,:,slice_idx));
colormap(gray);       % Use grayscale colormap for display
colorbar;             % Optional: display a colorbar to show intensity scale
axis image off;       % Maintain aspect ratio and hide axes
title(sprintf('Slice %d at TE1', slice_idx));


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

%% Visualizing as 3D Volume:

% Define the DICOM output directory
dicom_dir = 'DICOM_Output';
if ~exist(dicom_dir, 'dir')
    mkdir(dicom_dir);
end

% Combine the slices from TE1 and TE2 experiments
%all_slices = cat(3, img_TE1, img_TE2); % Concatenate along the 3rd dimension
% Use slices 1 to 14 (TE1)
%all_slices = cat(3, img_TE1);
selected_slices = img_TE1;
ReadoutFOV = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV;
PhaseFOV = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV;

NCol = twix.image.dataSize(1);
NLin = twix.image.dataSize(3);

% Define voxel dimensions (adjust if needed)
voxel_size_x = ReadoutFOV/NCol; % in mm
voxel_size_y = PhaseFOV/NLin; % in mm
voxel_size_z = sli_thickness; % slice thickness in mm

% Create a DICOM file for each slice
num_slices = size(selected_slices, 3);
for slice_idx = 1:num_slices
    % Extract the slice
    slice_data = selected_slices(:, :, slice_idx);
    
    % Create DICOM metadata
    metadata = struct();
    metadata.PatientName = 'PatientName'; % Replace with actual name if available
    metadata.PatientID = '12345'; % Replace with actual patient ID
    metadata.Modality = 'MR';
    metadata.Manufacturer = 'SIEMENS';
    metadata.InstitutionName = 'Tel-Aviv University';
    metadata.SliceThickness = voxel_size_z;
    metadata.PixelSpacing = [voxel_size_x, voxel_size_y];
    metadata.SliceLocation = slice_idx * voxel_size_z;
    metadata.SeriesDescription = 'Spin Echo MRI';
    metadata.StudyDescription = 'Brain MRI Study';
    % Add SOPClassUID to the metadata
    metadata.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4'; % MR Image Storage UID
    metadata.SOPInstanceUID = dicomuid(); % Generate a unique identifier for each slice
    metadata.SeriesInstanceUID = dicomuid(); % Unique identifier for the series
    metadata.StudyInstanceUID = dicomuid(); % Unique identifier for the study

    % Add InstanceNumber to keep slices ordered
    metadata.InstanceNumber = slice_idx;

    % Define the output DICOM file name
    dicom_filename = fullfile(dicom_dir, sprintf('Slice_%03d.dcm', slice_idx));
    
    % Write the DICOM file with the additional metadata
    dicomwrite(uint16(slice_data), dicom_filename, metadata, 'CreateMode', 'Copy');
end

% Visualize as a 3D volume
%[x, y, z] = meshgrid(1:size(selected_slices, 2), 1:size(selected_slices, 1), (1:size(selected_slices, 3)) * voxel_size_z);
%volumeViewer(selected_slices);

% Calculate the aspect ratio for rescaling
aspect_ratio = [voxel_size_y, voxel_size_x, voxel_size_z];

% Visualize as 3D volume with the correct aspect ratio
volumeViewer(selected_slices, 'ScaleFactors', aspect_ratio);

disp(['DICOM files saved to: ', dicom_dir]);