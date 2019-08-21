%% load parent folder %%

warning off

uiwait(msgbox('Load parent folder'));
parent_d = uigetdir('');

matlab_folder = cd;
cd(parent_d)
listing = dir('**/theta_nuclei_*.mat');
cd(matlab_folder)

%% open one file at a time and perform analysis %%

n_files = length(listing);
costheta_median = zeros(n_files,1);

for file_list = 1:n_files
    
    % file and directory name
    file = listing(file_list).name;
    directory = listing(file_list).folder;
    
    % output name
    slash_indeces = strfind(directory,'/');
    output_name = directory(slash_indeces(end)+1:end);
    
    % load theta nuclei orientation
    theta = load(fullfile(directory,file));
    theta = theta.theta;
    
    % calculate costheta similarity
    costheta_all = [];
    for ii = 1:length(theta)-1
        
        theta_ref = theta(ii,1);
        theta_temp = theta(ii+1:end,1);
        
        costheta = cosd(abs(theta_ref - theta_temp));
        costheta_all = [costheta_all; costheta];
        
        clear theta_temp costheta
        
    end
    
    % save all angles comparison
    save(fullfile(directory,['costheta_nuclei_', output_name,'.mat']), ...
        'costheta_all');
    
    % calculate median for each image
    costheta_median(file_list,1) = median(costheta_all);
    clear costheta_all
    
end

%% save %%
save(fullfile(parent_d, 'costheta_nuclei_median.mat'), ...
    'costheta_median')
