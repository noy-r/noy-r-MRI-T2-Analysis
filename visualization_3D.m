clc , clear, close all

%% Add the 'mapVBVD-main' folder to MATLAB's path %
addpath('mapVBVD-main');

% call mapVBVD directly 
twix = mapVBVD('meas_MID00172_FID189318_se.dat');
% In this tutorial we do not care for noise data, so let's get rid of the first measurement:
twix = twix{1};
twix.image.flagDoAverage = true; % automatically average the average dimension 
twix.image.flagRemoveOS = true; %% automatically removes oversampling

sli_thickness = twix.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness;

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

%% Visualizing as 3D Volume:

% Define the DICOM output directory
dicom_dir = 'DICOM_Output';
if ~exist(dicom_dir, 'dir')
    mkdir(dicom_dir);
end

% Combine the slices from TE1 and TE2 experiments
all_slices = cat(3, img_TE1, img_TE2); % Concatenate along the 3rd dimension
% Use slices 1 to 14 (TE1)
selected_slices = img_TE1;

% Define voxel dimensions (adjust if needed)
voxel_size_x = 1; % in mm
voxel_size_y = 1; % in mm
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
    metadata.SliceThickness = sli_thickness;
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

[x, y, z] = meshgrid(1:size(selected_slices, 2), 1:size(selected_slices, 1), (1:size(selected_slices, 3)) * voxel_size_z);
volumeViewer(selected_slices);

disp(['DICOM files saved to: ', dicom_dir]);