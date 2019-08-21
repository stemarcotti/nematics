%% input %%
d = uigetdir('');

slash_indeces = strfind(d,'/');
output_name = d(slash_indeces(end)+1:end);

file_name = ['orientation_matrix_', output_name,'.mat'];
orientation_matrix = load(fullfile(d, file_name));
orientation_matrix = orientation_matrix.Q; % [x, y, u, v]

% parameters
box_size_um = 50;	% [um]
mu2px = 1;          % [um]
box_size = round(box_size_um / mu2px);  % [px]

%% order parameter

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
        cos2theta = cos(2.*theta);
        sin2theta = sin(2.*theta);
        
        av_cos2theta = nanmean(cos2theta);
        av_sin2theta = nanmean(sin2theta);
        
        Qloc_temp = sqrt(av_cos2theta^2 + av_sin2theta^2);
        Qloc(k,:) = [x_ref, y_ref, Qloc_temp];
        
        clear cond1 cond2 cond3 cond4 cond
        clear data_temp_box
        clear u v theta
        clear cos2theta sin2theta
    end
    
end

Qloc(~any(Qloc,2),:) = [];  % delete zero rows

%% plot %%

spacing_matrix = 11;
size_matrix = 100;

Qloc1 = Qloc;

Qloc1(:,1) = Qloc1(:,1) ./ spacing_matrix;
Qloc1(:,2) = Qloc1(:,2) ./ spacing_matrix;

Qloc_matrix = zeros(size_matrix,size_matrix) .* NaN;

for ii = 1:size(Qloc1,1)
    
    idx = Qloc1(ii,2);
    idy = Qloc1(ii,1);
    Qloc_matrix(idx, idy) = Qloc1(ii,3);
    
end

figure
h = imshow(Qloc_matrix, []);
colormap('jet')
caxis([0, 1])
c = colorbar;

% black background
set(h, 'AlphaData', ~isnan(Qloc_matrix)) % set NaN to transparent
axis on;
set(gca, 'XColor', 'none', 'yColor', 'none', 'xtick', [], 'ytick', [], 'Color', 'k') % turn transparent to black
hold off
    
%% save %%

save(fullfile(d,['Qloc_', output_name,'.mat']), 'Qloc');
