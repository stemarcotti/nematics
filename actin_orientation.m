%% INPUT and parameters %%

[file_im, directory_im] = uigetfile('.tif');
slash_indeces = strfind(directory_im,'/');
output_name = directory_im(slash_indeces(end-1)+1:slash_indeces(end)-1);

window_size = 101;	% size of the sub windows (must be odd)
                    % use 21 for tilescans
                    % use 101 for actin stills

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

% mirror vectors
% this step is not needed in this case, as 
% atand(data(:,8) ./ data(:,7)) = atand(data(:,4) ./ data(:,3))
% kept this in for consistency with [fft_widowing.m]

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

figure
imshow(im)
hold on
quiver(data(:,1),data(:,2), data(:,7), data(:,8), 'g')

waitforbuttonpress
close

%% calculate vectors angle with horiziontal (-90 / +90)
theta = -atand(data(:,8) ./ data(:,7));
% negative sign here is to take care of difference with [regionprops]
% [nuclei_orientation.m] - angles in first quadrant (0 to +90) are
% considered positive, angles in fourth quadrant (0 to -90) are considered
% negative

save(fullfile(directory_im,['theta_actin_', output_name,'.mat']), 'theta');
clear