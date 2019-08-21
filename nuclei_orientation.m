%% %%
% load mask
[file,directory] = uigetfile('.tif');
nt = length(imfinfo([directory,'/', file]));

slash_indeces = strfind(directory,'/');
output_name = directory(slash_indeces(end-1)+1:slash_indeces(end)-1);
 
iptsetpref('ImshowBorder','tight');

%% %%
theta = [];
for k = 1:nt
    
    im = imread(fullfile(directory, file), k);
    im = im2double(im);
    
    bw = logical(im);
    bw_clean = bwareaopen(bw, 30);  % remove small objects (noise: area less than 30 px)
    bw_fill = imfill(bw_clean, 'holes');
    
    stats = regionprops(bw_fill, 'centroid', 'orientation');
    theta_temp = [stats(:).Orientation]';
    
    theta = [theta; theta_temp]; % nuclei major axis angle with horizontal [degrees]
end

figure
subplot(1,2,1)
imshow(bw_fill)
subplot(1,2,2)
polarhistogram(deg2rad(theta))

waitforbuttonpress
close 
%% %%
save(fullfile(directory,['theta_nuclei_', output_name,'.mat']), 'theta');
clear