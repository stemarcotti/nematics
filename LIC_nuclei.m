%% INPUT and parameters %%
uiwait(msgbox('Load actin image (.tif)'));
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
imshow(im)
hold on
quiver(data(:,1),data(:,2), data(:,3), data(:,4), 'g')

figure('Units', 'normalized', 'Position', [0.5 0.5 size(im,1) size(im,2)])
imshow(im)
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


%% %%

spacing_matrix = 51;
size_matrix = 19;

Q1 = Q;

Q1(:,1) = Q1(:,1) ./ spacing_matrix;
Q1(:,2) = Q1(:,2) ./ spacing_matrix;

Q_matrix_vx = zeros(size_matrix,size_matrix) .* NaN;
Q_matrix_vy = zeros(size_matrix,size_matrix) .* NaN;

for ii = 1:size(Q1,1)
    
    idx = Q1(ii,2);
    idy = Q1(ii,1);
    Q_matrix_vx(idx, idy) = Q1(ii,3);
    Q_matrix_vy(idx, idy) = Q1(ii,4);
    
end

%% interpolate %%

[X1,Y1] = meshgrid(linspace(1, size(Q_matrix_vx,1), size(im,1)));
Q_matrix_vx_int = interp2(Q_matrix_vx,X1,Y1);
Q_matrix_vy_int = interp2(Q_matrix_vy,X1,Y1);

cmap = ones(64, 3);

matlab_folder = cd;
cd('/Users/laste/Documents/OneDrive - King''s College London/git/LIC_NimaBigdelyShamlo')
previewvfield(Q_matrix_vy_int, Q_matrix_vx_int, cmap, size(im,1));
cd(matlab_folder)
hold on

%% get nuclei long axes %%

uiwait(msgbox('Load nuclei mask (.tif)'));
[file_im_nuclei, directory_im_nuclei] = uigetfile('.tif');

im_nuclei = imread(fullfile(directory_im_nuclei, file_im_nuclei));
im_nuclei = im2double(im_nuclei);

bw = logical(im_nuclei);
bw_clean = bwareaopen(bw, 30);  % remove small objects (noise: area less than 30 px)
bw_fill = imfill(bw_clean, 'holes');

stats = regionprops(bw_fill, 'centroid', 'orientation', 'majoraxislength');

x_centroid_n = zeros(length(stats),1);
y_centroid_n = zeros(length(stats),1);
for kk = 1:length(stats)
    x_centroid_n(kk,1) = stats(kk).Centroid(1);
    y_centroid_n(kk,1) = stats(kk).Centroid(2);
end
orientation_n = [stats(:).Orientation]';
axis_n = [stats(:).MajorAxisLength]';

x1_plot = x_centroid_n + (axis_n/2) .* cosd(-orientation_n);
x2_plot = x_centroid_n - (axis_n/2) .* cosd(-orientation_n);
y1_plot = y_centroid_n + (axis_n/2) .* sind(-orientation_n);
y2_plot = y_centroid_n - (axis_n/2) .* sind(-orientation_n);

for kk = 1:length(x1_plot)
    
    plot([x1_plot(kk,1),x2_plot(kk,1)], [y1_plot(kk,1),y2_plot(kk,1)], 'm-', 'linewidth', 3)
    
end

%% save %%
saveas(gcf, [directory_im '/LIC_nuclei_' output_name '.tif']);
clear