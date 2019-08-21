%% input %%
d = uigetdir('');

slash_indeces = strfind(d,'/');
output_name = d(slash_indeces(end)+1:end);

file_name = ['orientation_matrix_', output_name,'.mat'];
orientation_matrix = load(fullfile(d, file_name));
orientation_matrix = orientation_matrix.Q; % [x, y, u, v]

window_size = orientation_matrix(2,1) - orientation_matrix(1,1);
box_size = window_size * 2;

%% calculate order parameter %%

Qloc = zeros(size(orientation_matrix,1), 3);

for k = 1:size(orientation_matrix,1)
    
    % get ref vector data
    x_ref = orientation_matrix(k,1);
    y_ref = orientation_matrix(k,2);
    u_ref = orientation_matrix(k,3);
    v_ref = orientation_matrix(k,4);
    
    % remove border vectors from the analysis
    if ceil(x_ref - (box_size/2)) >= 0 && ...
            ceil(y_ref - (box_size/2)) >= 0 && ...
            floor(x_ref + (box_size/2)) <= max(orientation_matrix(:,1)) && ...
            floor(y_ref + (box_size/2)) <= max(orientation_matrix(:,2))
        
        % find data in box
        cond1 = orientation_matrix(:,1) <= floor(x_ref + box_size/2);
        cond2 = orientation_matrix(:,1) >= ceil(x_ref - box_size/2);
        cond3 = orientation_matrix(:,2) <= floor(y_ref + box_size/2);
        cond4 = orientation_matrix(:,2) >= ceil(y_ref - box_size/2);
        
        cond = cond1 & cond2 & cond3 & cond4;
        
        data_box = orientation_matrix(cond, :);
        
        % remove ref vector
        data_ref_idx = find(data_box(:,1) == x_ref & ...
            data_box(:,2) == y_ref);
        data_box(data_ref_idx, :) = [];
        
        % calculate angle between ref vector and vectors in box
        u = u_ref - data_box(:,3);
        v = v_ref - data_box(:,4);
        theta = atan(v./u);
        
        % calculate Qloc
        costheta2 = cos(theta).^2;
        av_costheta2 = nanmedian(costheta2);
        
        Qloc_temp = 2*(av_costheta2 - 0.5);
        Qloc(k,:) = [x_ref, y_ref, Qloc_temp];
        
        clear cond1 cond2 cond3 cond4 cond
        clear data_temp_box
        clear u v theta
        clear costheta2
    end
    
end

Qloc(~any(Qloc,2),:) = [];  % delete zero rows

%% make Qloc into matrix %%

size_matrix = max(Qloc(:,1)) / window_size;

Qloc1 = Qloc;

Qloc1(:,1) = Qloc1(:,1) ./ window_size;
Qloc1(:,2) = Qloc1(:,2) ./ window_size;

Qloc_matrix = zeros(size_matrix,size_matrix) .* NaN;

for ii = 1:size(Qloc1,1)
    
    idx = Qloc1(ii,2);
    idy = Qloc1(ii,1);
    Qloc_matrix(idx, idy) = Qloc1(ii,3);
    
end

%% interpolate %%

x = unique(Qloc(:,1))';
y = unique(Qloc(:,2))';
[X,Y] = meshgrid(x,y);

% set the points on which to interpolate
x_interp = min(x):max(x);
y_interp = min(y):max(y);
[X_interp,Y_interp] = meshgrid(x_interp,y_interp);

% interpolate
Qloc_matrix_interp = interp2(X,Y,Qloc_matrix,X_interp,Y_interp);

% smooth
Qloc_matrix_interp_smooth = imgaussfilt(Qloc_matrix_interp, 50);
imshow(Qloc_matrix_interp_smooth)
colormap('jet');

% save
saveas(gcf, [d '/heatmap_new_' output_name '.tif']);
