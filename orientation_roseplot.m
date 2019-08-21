%% load parent folder %%

warning off

uiwait(msgbox('Load parent folder'));
parent_d = uigetdir('');

matlab_folder = cd;
cd(parent_d)
listing = dir('**/theta_actin*.mat');
cd(matlab_folder)


%% open one file at a time and perform analysis %%

n_files = length(listing);

subplot_rows = 2;
subplot_n = length(listing);
range_polar = ([0:30:360] - 45) * pi / 180;

for file_list = 1:n_files
    
    % file and directory name
    file = listing(file_list).name;
    directory = listing(file_list).folder;
    
    % output name and cell ID
    slash_indeces = strfind(directory,'/');
    output_name = directory(slash_indeces(end)+1:end);
    
    theta_actin = load(fullfile(directory, file));
    theta_actin = theta_actin.theta;
    
    theta_nuclei = load(fullfile(directory, ['theta_nuclei_' output_name '.mat']));
    theta_nuclei = theta_nuclei.theta;
    
    subplot(subplot_rows, subplot_n / subplot_rows, file_list)
    polarhistogram(deg2rad(theta_nuclei), range_polar, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 3, ...
        'EdgeColor', 'k')
    hold on
    polarhistogram(deg2rad(theta_actin), range_polar, ...
        'Normalization', 'probability', ...
        'FaceColor', [0.9290 0.6940 0.1250])
   
end