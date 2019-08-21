%% load parent folder %%

warning off

uiwait(msgbox('Load parent folder'));
parent_d = uigetdir('');

matlab_folder = cd;
cd(parent_d)
listing = dir('**/*_mask.tif');
cd(matlab_folder)

%% open one file at a time and perform analysis %%

n_files = length(listing);
cell_shape_index = [];

for file_list = 1:n_files
    
    % file and directory name
    file = listing(file_list).name;
    directory = listing(file_list).folder;
    
    im = imread(fullfile(directory, file));
    im = im2double(im);
    
    bw = logical(im);
    bw_clean = bwareaopen(bw, 30);  % remove small objects (noise: area less than 30 px)
    bw_fill = imfill(bw_clean, 'holes');
    
    stats = regionprops(bw_fill, 'area', 'perimeter');
    
    % cell shape index [Versaevel et al. 2012, Nature Comm, Tiryaki et al. 2015, Cytometry]
    % CSI = (4?*area) / (perimeter^2)
    
    area_n = [stats(:).Area]';
    perimeter_n = [stats(:).Perimeter]';
    
    CSI = (4*pi * area_n) ./ (perimeter_n .^2);
    
    cell_shape_index = [cell_shape_index; CSI];
    
end

%% save %%

save(fullfile(parent_d,'n2_d1_shape_index.mat'), 'cell_shape_index');
