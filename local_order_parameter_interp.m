%% input %%
% matrix a is Fiji OrientationJ output

% assign variables
frame = a(:,3);
data = [a(:,1), a(:,2), a(:,4), a(:,5)]; % [x, y, u, v]

% parameters
box_size_um = 8;	% [um]
mu2px = 0.6;        % [um]
box_size = round(box_size_um / mu2px);  % [px]

%% order parameter %%

current_frame = 0;
while current_frame <= max(frame)
    
    % look at one frame at a time
    data_temp_idx = find(frame == current_frame);
    data_temp = data(data_temp_idx, :);
    
    %% interpolate
    % get meshgrid arrangement of points where vectors are calculated in Fiji
    x = unique(data_temp(:,1))';
    y = unique(data_temp(:,2))';
    [X,Y] = meshgrid(x,y);
    
    % order the column vectors [dx] and [dy] into matrix form
    for ii = 1:length(data_temp)
        dx(data_temp(ii,1), data_temp(ii,2)) = data_temp(ii,3);
    end
    dx(~any(dx,2),:) = [];  % delete zero rows
    dx(:,~any(dx,1)) = [];  % delete zero columns
    dx = dx';
    
    for ii = 1:length(data_temp)
        dy(data_temp(ii,1), data_temp(ii,2)) = data_temp(ii,4);
    end
    dy(~any(dy,2),:) = [];  % delete zero rows
    dy(:,~any(dy,1)) = [];  % delete zero columns
    dy = dy';
    
    % set the points on which to interpolate
    x_interp = min(x):max(x);
    y_interp = min(y):max(y);
    [X_interp,Y_interp] = meshgrid(x_interp,y_interp);
    
    % interpolate
    dx_interp = interp2(X,Y,dx,X_interp,Y_interp);
    dy_interp = interp2(X,Y,dy,X_interp,Y_interp);
    
    %% order parameter
    
    % correct shift
    X_interp_shift = X_interp - min(x_interp);
    Y_interp_shift = Y_interp - min(y_interp);

    % create input matrix
    X_interp_vect = X_interp_shift(:);
    Y_interp_vect = Y_interp_shift(:);
    dx_interp_vect = dx_interp(:);
    dy_interp_vect = dy_interp(:);
    data_temp_interp = [X_interp_vect, Y_interp_vect, dx_interp_vect, dy_interp_vect];
    
    
    Qloc = zeros(max(data_temp_interp(:,1)), max(data_temp_interp(:,2)));
    for k = 1:size(data_temp_interp,1)
        
        % get ref vector data
        x_ref = data_temp_interp(k,1);
        y_ref = data_temp_interp(k,2);
        u_ref = data_temp_interp(k,3);
        v_ref = data_temp_interp(k,4);
        
        % remove border vectors from the analysis
        if ceil(x_ref - (box_size/2)) >= 0 && ...
                ceil(y_ref - (box_size/2)) >= 0 && ...
                floor(x_ref + (box_size/2)) <= max(data_temp_interp(:,1)) && ...
                floor(y_ref + (box_size/2)) <= max(data_temp_interp(:,2))
            
            % find data in box
            cond1 = data_temp_interp(:,1) <= floor(x_ref + box_size/2);
            cond2 = data_temp_interp(:,1) >= ceil(x_ref - box_size/2);
            cond3 = data_temp_interp(:,2) <= floor(y_ref + box_size/2);
            cond4 = data_temp_interp(:,2) >= ceil(y_ref - box_size/2);
            
            cond = cond1 & cond2 & cond3 & cond4;
            
            data_temp_box = data_temp_interp(cond, :);
            
            % remove ref vector
            data_ref_idx = find(data_temp_box(:,1) == x_ref & ...
                data_temp_box(:,2) == y_ref);
            data_temp_box(data_ref_idx, :) = [];
            
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
    
    
    
end







