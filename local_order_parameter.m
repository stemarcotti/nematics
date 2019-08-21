%% input
% matrix a is Fiji OrientationJ output

% assign variables
frame = a(:,3);
data = [a(:,1), a(:,2), a(:,4), a(:,5)]; % [x, y, u, v]

% parameters
box_size_um = 50;	% [um]
mu2px = 1;        % [um]
box_size = round(box_size_um / mu2px);  % [px]

%% order parameter
current_frame = 0;
while current_frame <= max(frame)
    
    % look at one frame at a time
    data_temp_idx = find(frame == current_frame);
    data_temp = data(data_temp_idx, :);
    
    Qloc = zeros(max(data_temp(:,1)), max(data_temp(:,2)));
    for k = 1:size(data_temp,1)
        
        % get ref vector data
        x_ref = data_temp(k,1);
        y_ref = data_temp(k,2);
        u_ref = data_temp(k,3);
        v_ref = data_temp(k,4);
        
        % remove border vectors from the analysis
        if ceil(x_ref - (box_size/2)) >= 0 && ...
                ceil(y_ref - (box_size/2)) >= 0 && ...
                floor(x_ref + (box_size/2)) <= max(data_temp(:,1)) && ...
                floor(y_ref + (box_size/2)) <= max(data_temp(:,2))
            
            % find data in box
            cond1 = data_temp(:,1) <= floor(x_ref + box_size/2);
            cond2 = data_temp(:,1) >= ceil(x_ref - box_size/2);
            cond3 = data_temp(:,2) <= floor(y_ref + box_size/2);
            cond4 = data_temp(:,2) >= ceil(y_ref - box_size/2);
            
            cond = cond1 & cond2 & cond3 & cond4;
            
            data_temp_box = data_temp(cond, :);
            
            % remove ref vector
            data_ref_idx = find(data_temp_box(:,1) == x_ref & ...
                data_temp_box(:,2) == y_ref);
            data_temp_box(data_ref_idx, :) = [];
            
            %% CHECK THIS BIT %% 
            % calculate angle between ref vector and vectors in box
            u = u_ref - data_temp_box(:,3);
            v = v_ref - data_temp_box(:,4);
            theta = atan(v./u);
            
            % calculate Qloc
            cos2theta = cos(2.*theta);
            sin2theta = sin(2.*theta);
            
            av_cos2theta = nanmean(cos2theta);
            av_sin2theta = nanmean(sin2theta);
            
            Qloc(x_ref, y_ref) = sqrt(av_cos2theta^2 + av_sin2theta^2);
            
            clear cond1 cond2 cond3 cond4 cond
            clear data_temp_box
            clear u v theta
            clear cos2theta sin2theta
        end
        
    end
    
    
    current_frame = current_frame+1;
end