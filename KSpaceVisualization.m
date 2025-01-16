clc , clear, close all

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

