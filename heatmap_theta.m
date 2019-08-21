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

%% angles heatmap

theta = [data(:,1), data(:,2), data(:,5)];
theta_rotated = [data(:,1), data(:,2), data(:,6)];

interval_size = theta(2,1) - theta(1,1);
size_matrix = max(theta(:,1)) / interval_size;

T = theta;
Tr = theta_rotated;

T(:,1) = T(:,1) ./ interval_size;
T(:,2) = T(:,2) ./ interval_size;

Tr(:,1) = Tr(:,1) ./ interval_size;
Tr(:,2) = Tr(:,2) ./ interval_size;

T_matrix = zeros(size_matrix,size_matrix) .* NaN;
Tr_matrix = zeros(size_matrix,size_matrix) .* NaN;

for ii = 1:size(T,1)
    
    idx = T(ii,2);
    idy = T(ii,1);
    T_matrix(idx, idy) = T(ii,3);
    
end

for ii = 1:size(Tr,1)
    
    idx = Tr(ii,2);
    idy = Tr(ii,1);
    Tr_matrix(idx, idy) = Tr(ii,3);
    
end

%% interpolate %%

x = unique(theta_rotated(:,1))';
y = unique(theta_rotated(:,2))';
[X,Y] = meshgrid(x,y);

% set the points on which to interpolate
x_interp = min(x):max(x);
y_interp = min(y):max(y);
[X_interp,Y_interp] = meshgrid(x_interp,y_interp);

% interpolate
Tr_matrix_interp = interp2(X,Y,Tr_matrix,X_interp,Y_interp);
T_matrix_interp = interp2(X,Y,T_matrix,X_interp,Y_interp);

% smooth
Tr_matrix_interp_smooth = imgaussfilt(Tr_matrix_interp, 20);
T_matrix_interp_smooth = imgaussfilt(T_matrix_interp, 20);

% figure
% imshow(T_matrix_interp_smooth, [])
% colormap('hsv');
% caxis([-90, 90])
% saveas(gcf, [directory_im '/heatmap_theta_-9090_' output_name '.tif']);
% close



figure

a = colormap('hsv');
b = flipud(a);
a_length = length(a);
c = a;

c(a_length/2+1:end, :) = b(a_length/2+1:end, :);

% a = colormap('hsv');
% a = a(1:2:end, :);
% b = flipud(a);
% c = [a;b];

imshow(Tr_matrix, [])
colormap(c);
caxis([0, 180])
saveas(gcf, [directory_im '/heatmap_theta_0180_wrapped_nointerp_' output_name '.tif']);


clear; clc