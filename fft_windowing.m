%% INPUT and parameters %%

[file_im, directory_im] = uigetfile('.tif');
slash_indeces = strfind(directory_im,'/');
output_name = directory_im(slash_indeces(end-1)+1:slash_indeces(end)-1);

window_size = 55;	% size of the sub windows (must be odd)
                    % use 21 for tilescans
                    % use 101 for actin stills
                    % use 55 for heatmaps

overlap = 0.5;      % value should be between 0 (complete overlap) and 1 (no overlap)
st = 2 ;            % window of size 2*st+1 to compare vectors to (order parameter calculations)
checkpoint = 0;     % threshold sum in each window to do the calculation
mask_method = 1;	% global = 1; local = 2
figures = 0;        % plot figures after calculation
maskname = [];      % analyzes the entire image if empty

% load image
im = double(imread(fullfile(directory_im, file_im))) / 255;

%% calculate orientation matrix [FFTAlignment.m from Cetera et al. Nat Commun 2014]

vectors_image_filename = [directory_im, 'vectors_', file_im];
vectors_histogram_image_filename = [directory_im, 'vectors_histogram_', file_im];
FFT_alignment_data = FFTAlignment(im, window_size, overlap, st, checkpoint,...
    mask_method, maskname, figures, vectors_image_filename,...
    vectors_histogram_image_filename);

data = [FFT_alignment_data.pos, FFT_alignment_data.vec];

%% mirror vectors

data(:,5) = atan2d(data(:,4), data(:,3));
for k = 1:size(data,1)
    if data(k,5) >= 0
        data(k,6) = data(k,5);
    else
        data(k,6) = 180+data(k,5);
    end
end
data(:,7) = cosd(data(:,6));
data(:,8) = sind(data(:,6));

figure('Units', 'normalized', 'Position', [0.25 0.5 size(im,1) size(im,2)])
imshow(im, [])
hold on
quiver(data(:,1),data(:,2), data(:,3), data(:,4), 'g')

figure('Units', 'normalized', 'Position', [0.5 0.5 size(im,1) size(im,2)])
imshow(im, [])
hold on
quiver(data(:,1),data(:,2), data(:,7), data(:,8), 'm')

waitforbuttonpress

mirror_q = questdlg('Do you want to mirror the vectors', ...
    'Mirror?',...
    'No (left)', 'Yes (right)', 'No (left)');

close all

if strcmp(mirror_q, 'Yes (right)')
    Q = [data(:,1),data(:,2), data(:,7), data(:,8)];
else
    Q = [data(:,1),data(:,2), data(:,3), data(:,4)];
end

%% plot and save

save(fullfile(directory_im, ['orientation_matrix_', output_name,'.mat']), ...
    'Q');

clear; clc; close all
