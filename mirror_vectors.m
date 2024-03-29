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

%% call FFTAlignment to obtain vector field [-90�;90�] %%
vectors_image_filename = [directory_im, 'vectors_', file_im];
vectors_histogram_image_filename = [directory_im, 'vectors_histogram_', file_im];
FFT_alignment_data = FFTAlignment(im, window_size, overlap, st, checkpoint,...
    mask_method, maskname, figures, vectors_image_filename,...
    vectors_histogram_image_filename);

% data [position_x position_y
%       u_-90/90 v_-90/90]
data = [FFT_alignment_data.pos, FFT_alignment_data.vec];

%% obtain theta from vector field and save vector components in [0�;180�] %%

% data [position_x position_y
%       u_-90/90 v_-90/90
%       theta_-90/90 theta_0-180
%       u_0-180(unit) v_-180(unit)]

data(:,5) = atan2d(data(:,4), data(:,3));   % theta [-90�;90�]
for k = 1:size(data,1)
    if data(k,5) >= 0
        data(k,6) = data(k,5);      % theta [0�;180�]
    else
        data(k,6) = 180+data(k,5);  % theta [0�;180�]
    end
end
data(:,7) = cosd(data(:,6));    % u_0-180(unit vector)
data(:,8) = sind(data(:,6));    % v_0-180(unit vector)

%% perform selective vector rotation %%

% A: original vector field [-90�;90�]	-> I and IV quadrants
% B: rotated vector field [0�;180�]     -> I and II quadrants

% rationale
% * vectors in the I quadrant have no problems (as no discontinuity neither in A or B)
% * A as problems with vertical switch, B with horizontal switch (boundary discontinuity)
% * decide whether rotating each vector in II and IV quadrant depending on theta of closest neighbour(s) in I quadrant

% number of neighbours to take into account
neighbours = 5;

% find vectors in I quadrant (these are the same in A and B)
idx_AB = find(data(:,5) == data(:,6));
AB = data(idx_AB, :);

% find vectors in II and IV quadrants (these are NOT the same in A and B)
idx_notAB = find(data(:,5) ~= data(:,6));
notAB = data(idx_notAB, :);

% initialise matrix to store rotated vectors 
notAB_new = zeros(size(notAB,1),4);

% for each vector in II and IV quadrants
for k = 1:size(notAB,1)
    
    % [temp] is current vector to be rotated
    temp = notAB(k,:);
    temp_pos = temp(:, 1:2);
    
    % calculate distance of current vectors to all vectors in I quadrant
    dist = sqrt((temp_pos(1) - AB(:,1)).^2 + (temp_pos(2) - AB(:,2)).^2);
    % find closest n neighbours
    [~, min_idx] = mink(dist, neighbours);
    
    % average theta of neighbours to decide if rotating
    if mean(AB(min_idx, 5)) <= 45
        temp_vect = notAB(k,3:4) ./ norm(notAB(k,3:4)); % unit vector
        notAB_new(k,:) = [notAB(k,1:2) temp_vect];      % do not rotate if average theta of neighbours is between 0� and 45� (close to IV quadrant) 
    else
        notAB_new(k,:) = notAB(k,[1,2,7,8]);            % rotate if average theta of neighbours is between 45� and 90� (close to II quadrant) 
    end
    
end

% save unit vectors also for I quadrant vectors
AB_new = AB(:,[1,2,7,8]);

% [data_new] contains I quadrant vectors + selectively rotated ones
data_new = [AB_new; notAB_new];


%% PLOT %%

% original vector field [-90�;90�]
figure
imshow(im)
hold on
quiver(data(:,1),data(:,2), data(:,3), data(:,4), 'g', 'linewidth', 1)
pause(1)

% rotated vector field [0�;180�]
figure
imshow(im)
hold on
quiver(data(:,1),data(:,2), data(:,7), data(:,8), 'm', 'linewidth', 1)
pause(1)

% original and rotated together
figure
imshow(im)
hold on
quiver(data(:,1),data(:,2), data(:,3), data(:,4), 'g', 'linewidth', 1)
quiver(data(:,1),data(:,2), data(:,7), data(:,8), 'm', 'linewidth', 1)
pause(1)

% selectively rotated vector field
figure
imshow(im)
hold on
quiver(data_new(:,1),data_new(:,2), data_new(:,3), data_new(:,4), 'w', 'linewidth', 1)
pause(1)