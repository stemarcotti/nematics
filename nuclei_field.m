%% %%
% load mask
[file,directory] = uigetfile('.tif');
nt = length(imfinfo([directory,'/', file]));

%% %%
for k = 1:nt
    
    im = imread(fullfile(directory, file), k);
    im = im2double(im);
    
    bw = logical(im);
    bw_clean = bwareaopen(bw, 2);
    
    bw_dilate = imdilate(bw_clean, strel('disk', 3));
    bw_dilate_clean = bwareaopen(bw_dilate, 300);
    bw_fill = imfill(bw_dilate_clean, 'holes');
    bw_erode = imerode(bw_fill, strel('disk', 7));
    
    [L,~] = bwlabel(bw_erode);
    stats = regionprops(L, 'Centroid', 'Orientation', 'Area');
    
    data_centroid = cat(1, stats.Centroid);
    data_area = cat(1, stats.Area);
    data_orientation = cat(1, stats.Orientation);
    
    %     data = [data_centroid, data_area, data_orientation];
    %
    %     q3 = quantile(data_area,0.75);
    %     q1 = quantile(data_area,0.25);
    %     whisker = q3 + 1.5 * (q3 - q1);
    %     idx = find(data_area > whisker);
    %     data(idx,:) = [];
    
    
    imshow(bw_erode)
    hold on
    x = data_centroid(:,1);
    y = data_centroid(:,2);
    
    %     u = zeros(size(x,1),1);
    %     v = zeros(size(x,1),1);
    %     for kk = 1:size(x,1)
    %         alpha = data_orientation(kk,1);
    %         if (alpha >= 0 && alpha < 30) || (alpha < 0 && alpha >= -30)
    %             u(kk,1) = data_box(kk,1)+data_box(kk,3)-data_centroid(kk,1);
    %             v(kk,1) = 0;
    %         elseif (alpha >= 30 && alpha < 60)
    %             u(kk,1) = data_box(kk,1)+data_box(kk,3)-data_centroid(kk,1);
    %             v(kk,1) = data_box(kk,2)-data_centroid(kk,2);
    %         elseif (alpha >= 60) || (alpha < -60)
    %             u(kk,1) = 0;
    %             v(kk,1) = data_box(kk,2)-data_centroid(kk,2);
    %         else
    %             u(kk,1) = data_box(kk,1)-data_centroid(kk,1);
    %             v(kk,1) = data_box(kk,2)-data_centroid(kk,2);
    %         end
    %
    %     end
    
    u = cosd(data_orientation);
    v = sind(data_orientation);
    for kk = 1:size(x,1)
        alpha = data_orientation(kk,1);
        if alpha < -30
            u(kk,1) = -u(kk,1);
            v(kk,1) = -v(kk,1);
        end
    end

    quiver(data_centroid(:,1), data_centroid(:,2), cosd(data_orientation), sind(data_orientation), 'g')
    quiver(x,y,u,v, 'm')
    
%     field_u = zeros(size(x,1),size(y,2));
%     field_v = zeros(size(x,1),size(y,2));
%     for kk = 1:size(x,1)
%         for jj = 1:size(y,1)
%             field_u(kk,jj) = u(kk,1);
%             field_v(kk,jj) = v(kk,1);
%         end
%     end
        
    data_temp = [round(x) round(y) u v];
    
    % order the column vectors [dx] and [dy] into matrix form
    for ii = 1:size(data_temp,1)
        dx(data_temp(ii,1), data_temp(ii,2)) = data_temp(ii,3);
    end
    dx(~any(dx,2),:) = [];  % delete zero rows
    dx(:,~any(dx,1)) = [];  % delete zero columns
    dx = dx';
    
    for ii = 1:size(data_temp,1)
        dy(data_temp(ii,1), data_temp(ii,2)) = data_temp(ii,4);
    end
    dy(~any(dy,2),:) = [];  % delete zero rows
    dy(:,~any(dy,1)) = [];  % delete zero columns
    dy = dy';
    
    % set the points on which to interpolate
    x_grid = 1:size(dx,1);
    y_grid = 1:size(dx,2);
    [X,Y] = meshgrid(x_grid,y_grid);
    
    x_interp = 1:size(im,1);
    y_interp = 1:size(im,2);
    [X_interp,Y_interp] = meshgrid(x_interp,y_interp);
    
    % interpolate
    dx_interp = interp2(X,Y,dx,X_interp,Y_interp);
    dy_interp = interp2(X,Y,dy,X_interp,Y_interp);
    
    
    

end