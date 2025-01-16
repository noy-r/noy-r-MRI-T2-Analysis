clc , clear, close all

%% Add the 'mapVBVD-main' folder to MATLAB's path %
addpath('mapVBVD-main');

% call mapVBVD directly 
twix = mapVBVD('meas_MID00172_FID189318_se.dat');
% In this tutorial we do not care for noise data, so let's get rid of the first measurement:
twix = twix{1};
twix.image.flagDoAverage = true; % automatically average the average dimension 
twix.image.flagRemoveOS = true; %% automatically removes oversampling

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
TEs = twix.hdr.MeasYaps.alTE;  % cell array containing TEs in microseconds
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