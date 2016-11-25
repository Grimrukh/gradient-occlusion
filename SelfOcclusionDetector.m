function [FINAL_DISPLAY] = SelfOcclusionDetector(image,mask)
    
    % Detects luminance falloff correlations along detected edges. Returns
    % local correlations (red) and the best fit of global correlation
    % (green).
    
    %% Initialize parameters
    
    profile on
    
    if nargin < 2
        mask = [];
    end
    
    SHOW_MASK = 1;
    
    Rsquare_threshold = 0.3;            % Minimum rSq for image display
    inner_Rsquare_threshold = 0.6;      % Minimum rSq for use in inner correlation calculation
    slope_threshold = 0.1;              % Excludes "perfect" correlations with no gradients
    
    probe_distances = 3:5;              % Ring distances from edge
    probe_directions = [-pi/2 pi/2];    % Ring directions from edge
    
    
    
    %% MAIN FUNCTION    
    
    [image,bit_depth] = ReadImage(image);
    image = AutoContrast(image,bit_depth);
    image = double(image);
    imsize = size(image);
    grey = mean(image(:));
    
    if nargin > 1
        [binary,~,vectors,edgelist] = EdgesFromMask(mask);
    else
        [~,vectors,edgelist] = EdgesFromImage(image);
    end
    
    edgelist = BreakEdges(edgelist,100);
    
    edge_numbers = 1:length(edgelist);
    array_size = [length(edge_numbers), length(probe_directions), length(probe_distances)];

    
    EDGE_DATA = GetShiftedPixels(edgelist,vectors);
    
    % CORRELATIONS
    
    local_image = ComputeLocalCorrelation(edgelist,EDGE_DATA);
    global_image = ComputeGlobalCorrelation(edgelist,EDGE_DATA);
    
    if SHOW_MASK
        size(image)
        image_RGB = repmat(image,[1 1 3]);
        mask_binary = repmat(binary,[1 1 3]);
        size(mask_binary)
        size(image_RGB)
        masked_image = (mask_binary*grey) + (~mask_binary.*double(image_RGB));
        local_image = (mask_binary*grey) + (~mask_binary.*double(local_image));
        global_image = (mask_binary.*grey) + (~mask_binary.*double(global_image));        
    end
    
    
    % FINAL IMAGE DISPLAY
    
    FINAL_DISPLAY = uint16([image_RGB masked_image; ...
                     local_image global_image]);
    figure(20), imshow(FINAL_DISPLAY);
    
    profile viewer
    
    
    
            
    %% Read intensity image
    
    function [im_out,bit_depth] = ReadImage(im)
        
        if isa(im,'char')
            try
                im = ReadGray(strcat('images/',im,'.tif'));
            catch
                error('Could not find or read image file.');
            end
        end
        
        if isa(im,'uint8')
            bit_depth = 8;
        elseif isa(im,'uint16')
            bit_depth = 16;
        else
            error('Image is not an integer type.');
        end
        
        if size(im,3) > 1
            im_out = rgb2gray(im);
        else
            im_out = im;
        end
        
    end
    
    
    %% Scale intensity image
    
    function im_out = AutoContrast(im,bits)
       
        im = double(im);
        im_max = max(im(:));
        im_min = min(im(:));
        
        if bits == 8
            im_out = uint8( 2^bits*(im-im_min) / (im_max-im_min) );
        elseif bits == 16
            im_out = uint16( 2^bits*(im-im_min) / (im_max-im_min) );
        else
            error('Invalid bit-depth for scaling.');
        end
        
    end
    
    
    %% Read mask image
    
    function [mask_out] = ReadMask(mask_im)
        
        if isa(mask_im,'char')
            try
                mask_im = ReadGray(strcat('images/masks/',mask_im,'.tif'));
            catch
                error('Could not find or read image file.');
            end
        end  
        
        mask_out = imbinarize(mask_im); 
        
    end
    

   
    %% Create edge list from an input mask image
    
    function [m,edgs,vecs,edgelist] = EdgesFromMask(mask_image)        
        
        m = ReadMask(mask_image);
        [edgs,vecs] = AnalyseMaskEdges(m,[]);
        edgelist = edgelink(edgs);  
        
    end
    
    
    
    %% Detect edges
    
    function edges = DetectEdges(image,threshold)
        
        if nargin<2
            threshold = [];
        end
        
        edges = edge(image,'Canny',threshold);
        dilated = imdilate(edges,strel('disk,3'));
        edges = bwmorph(dilated,'thin',Inf);
        
    end
    
    %% Create edge lists
    
    function [edges,vectorimages,edgelist] = EdgesFromImage(image,threshold)
    
        lowpass = imgaussfilt(image,5);
        highpass = image - lowpass;

        if length(threshold) > 1            
            % Try all Canny thresholds in list and select the one with the 
            % highest average edge length after tiny edges are removed.
            mean_edge_lengths = zeros(length(threshold));
            edgelists(length(threshold)).edgelist = [];

            for t = 1:length(threshold)
                edges = DetectEdges(highpass,threshold(t));                
                [edgelist, ~] = edgelink(edges);                
                edgelists(t).edgelist = edgelist;
                edgelists(t).edge_image = edges;
                edge_lengths = cellfun('length',edgelist);
                mean_edge_lengths(t) = nanmean(edge_lengths);
            end

            [~,max_mean_edge_length_index] = max(mean_edge_lengths);
            edgelist = edgelists(max_mean_edge_length_index).edgelist;

        else            
            edges = DetectEdges(highpass);
            [edgelist, ~] = edgelink(edges);            
        end
        
        vectorimages(length(edgelist)) = [];
        
        for e = 1:length(edgelist)                        
            if size(edgelist{e},1) > 10
                vectorimages{e} = deg2rad(skeletonOrientation(edges,9,1,edgelist{e}));
            end            
        end
        
    end
    
    %% Break up edges over limit
    
    function edgelist = BreakEdges(edgelist,max_size)
    
        % Calculate size of new edgelist.
    
        edge_lengths = cellfun('length',edgelist);
        current_edgelist_length = length(edgelist);
        new_edgelist_length = sum(ceil(edge_lengths/100));
        edgelist{1,new_edgelist_length} = [];
        
        % Break up edgelist.        
        
        write_index = current_edgelist_length + 1;
        for e = 1:current_edgelist_length

            edge_subs = edgelist{e};

            while size(edge_subs,1) > max_size

                edgelist{write_index} = edge_subs(1:max_size,:);
                edgelist{e} = edge_subs(max_size+1:end,:);
                edge_subs = edgelist{e};
                write_index = write_index + 1;

            end
        end
    end
        
    function EDGE_DATA = GetShiftedPixels(edgelist,vectors)
    
        % Initialize edge data        
        EDGE_DATA(length(edgelist)).Intensity = 0;

        %% Calculate orientation of each edge and fetch shifted pixels

        for e = edge_numbers

            current_edge_subs = edgelist{e};
            edge_length = size(current_edge_subs,1);
            
            % Skip very short edges (no vectors for them)
            if edge_length <= 10
                continue
            end
            
            % Reset data fields
            current_edge_intensity = zeros(edge_length,length(probe_distances),length(probe_directions));
            current_edge_angle = zeros(edge_length,1);
            current_edge_shifted = zeros(edge_length,length(probe_distances),length(probe_directions));
            
            for pix = 1:edge_length

                pix_row = current_edge_subs(pix,1);
                pix_col = current_edge_subs(pix,2);

                % Reset shifted pixel fields
                pix_shifted = zeros(length(probe_distances),length(probe_directions));
                pix_shifted_int = pix_shifted;

                % Create shift vectors
                probe_angles = vectors(pix_row,pix_col) + probe_directions;
                shift_col = round(probe_distances' * cos(probe_angles));
                shift_row = round(probe_distances' * -sin(probe_angles));

                if (any(pix_row+shift_row(:) <= imsize(1)) && any(pix_col+shift_col(:) <= imsize(2)))

                    % Try grabbing each pixel.
                    for i = 1:length(probe_distances)
                        for j = 1:length(probe_directions)
                            try
                                pix_shifted(i,j) = sub2ind(imsize,pix_row+shift_row(i,j),pix_col+shift_col(i,j));
                                pix_shifted_int(i,j) = image(pix_shifted(i,j));
                            catch
                                % Do nothing.
                            end
                        end
                    end

                end

                current_edge_intensity(pix,:,:) = pix_shifted_int;
                current_edge_angle(pix) = vectors(pix_row,pix_col);
                current_edge_shifted(pix,:,:) = pix_shifted;

            end

            EDGE_DATA(e).Intensity = current_edge_intensity;
            EDGE_DATA(e).Orientation = current_edge_angle;
            EDGE_DATA(e).ShiftedPoints = current_edge_shifted;

        end
       
    end
    
    % All pixels of all edges have been checked. 
    
    % First, each edge's local correlation will be checked. The first
    % display image will show edge correlation strength based on individual
    % falloffs (with individual max angles).
    
    % Second, the average significant max_angle will enforce itself 
    % upon all the other edges and directions. This will give a sense of 
    % the global correlation.
    
    function [local_image,inner_image,angle_image] = ComputeLocalCorrelation(edgelist,EDGE_DATA)
    
        RUNNING_MAX_ANGLES = NaN * zeros(length(edgelist),1);
        RUNNING_MAX_INTENSITY = NaN * zeros(length(edgelist),1);
        RUNNING_RSQUARE = NaN * zeros(length(edgelist),1);
        RUNNING_EDGE_LENGTH = NaN * zeros(length(edgelist),1);
        color_image = zeros(imsize);
        angle_image = repmat(image,[1 1 3]);
        array_shift = imsize(1)*imsize(2);
        inner_color_image = zeros(imsize);
        
        disp('Running local edgewise correlations...');

        for e = edge_numbers

            % Perform separate correlation for every 'ring' (combination of
            % direction and distance).

            if size(edgelist{e},1) <= 10
                continue % Skip small edges
            end

            for curr_dir = 1:length(probe_directions)

                regression_slopes = zeros(length(probe_distances),1);
                regression_Rsquares = zeros(length(probe_distances),1);

                for curr_dst = 1:length(probe_distances)

                    % Reset Rsquare and B
                    slope = 0;
                    Rsquare = 0;

                    current_points = EDGE_DATA(e).ShiftedPoints(:,curr_dst,curr_dir);
                    current_intensity = EDGE_DATA(e).Intensity(:,curr_dst,curr_dir);

                    % Find maximum intensity of edge, and get corresponding angle.
                    [max_int, max_int_ind] = max(current_intensity);
                    [~, min_int_ind] = min(current_intensity);
                    max_int_angle = EDGE_DATA(e).Orientation(max_int_ind);
                    min_int_angle = EDGE_DATA(e).Orientation(min_int_ind);

                    current_intensity = current_intensity(current_points>0);
                    current_angles = EDGE_DATA(e).Orientation(current_points>0);
                    current_points = current_points(current_points>0);

                    new_color = 0;
                    angle_color = [0 0 0];

                    for ANG = [max_int_angle min_int_angle]

                        current_ang_dst = abs(ANG - current_angles);
                        pi_check = current_ang_dst > pi;

                        current_ang_dst = ~pi_check .* current_ang_dst + ...
                            pi_check .* (2*pi - current_ang_dst);

                        [O,ind] = sort(current_ang_dst);
                        G = current_intensity(ind);
                        P = current_points(ind);

                        if length(P) >= 10

                            X = [O ones(size(P))];
                            [~,R,~] = qr(X,0);
                            p = sum(abs(diag(R)) > max(length(P),2)*eps(R(1)));

                            if p >= 2                            
                                [B,~,~,~,stats] = regress(G,X);
                            else
                                B = 0;
                                stats = 0;
                            end

                            % Continue if this is the new best Rsquare and the
                            % slope is above threshold.
                            if (stats(1) > Rsquare && abs(B(1)) > slope_threshold)
                                Rsquare = stats(1);
                                if ANG == max_int_angle
                                    slope = B(1);
                                elseif ANG == min_int_angle
                                    slope = -B(1);
                                end

                                % Update images if Rsquare is above threshold.
                                if Rsquare > Rsquare_threshold
                                    new_color = (Rsquare - Rsquare_threshold) / (1 - Rsquare_threshold);
                                    angle_color = 255 * ang2rgb(ANG,0.5);
                                    if curr_dst == 1
                                        % Add to running max angles if this is
                                        % the first distance.
                                        if ANG == max_int_angle
                                            RUNNING_MAX_ANGLES(e) = max_int_angle;
                                            RUNNING_MAX_INTENSITY(e) = max_int;
                                            RUNNING_RSQUARE(e) = Rsquare;
                                            RUNNING_EDGE_LENGTH(e) = length(P);
                                        elseif ANG == min_int_angle
                                            % Minimize contribution to max angle.
                                            RUNNING_MAX_ANGLES(e) = min_int_angle + pi;
                                            RUNNING_MAX_INTENSITY(e) = 0.1;
                                            RUNNING_RSQUARE(e) = 0.1;
                                            RUNNING_EDGE_LENGTH(e) = 0.1;
                                        end
                                    end
                                else
                                    angle_color = [0 0 0];
                                    new_color = 0;
                                end                        
                            end
                        else
                            Rsquare = 0;
                            slope = 0;
                            new_color = 0;
                            angle_color = [0 0 0];
                        end
                    end

                    regression_slopes(curr_dst) = slope;
                    regression_Rsquares(curr_dst) = Rsquare;

                    angle_image(current_points) = angle_color(1);
                    angle_image(current_points+array_shift) = angle_color(2);
                    angle_image(current_points+2*array_shift) = angle_color(3);
                    color_image(current_points) = max(color_image(current_points),new_color);
                end

                % Run inner regression for curvature analysis

                if (length(P) >= 10 && mean(regression_Rsquares) > inner_Rsquare_threshold)
                    [BI,~,~,~,statsI] = regress(regression_slopes,[probe_distances' ones(size(probe_distances'))]);
                    RsquareI = statsI(1);                
                else
                    RsquareI = 0;
                    BI = 0;
                end

                % Set image color weight using Rsquare if the slope is
                % positive.

                if BI(1) > 0
                    new_color = RsquareI;
                else
                    new_color = 0;
                end

                inner_color_image(current_points) = max(inner_color_image(current_points),new_color);

            end

        end

        % Put color image on original image

        angle_image = uint8(angle_image);
        local_image = RgbBlend(image, color_image, 'red');       
        inner_image = RgbBlend(image, inner_color_image, 'blue');
        
    end
    
    function [global_image] = ComputeGlobalCorrelation(edgelist,EDGE_DATA)
    
        %% New global method: try all eight angles (sign ambiguous) for best global correlation.

        disp('Running global edgewise correlations...');

        max_angles = 0:pi/8:7*pi/8+0.01;
        
        % Keep track of properties of best rSquare

        RUNNING_max_avg_rSq = 0;
        RUNNING_max_rSq = NaN * zeros(array_size);
        RUNNING_ang = 0;
        RUNNING_pixels = [];

        for ang = max_angles

            % Reset Rsquares and pixel locations.
            
            [rSq,pixels,ring_lengths] = AnalyseGlobal(edgelist,EDGE_DATA,ang);

            % Weighted average of ALL current Rsquares using ring lengths.

            max_avg_rSq = wmean(rSq,ring_lengths,0);

            if max_avg_rSq > RUNNING_max_avg_rSq
                RUNNING_max_avg_rSq = max_avg_rSq;
                RUNNING_max_rSq = rSq;
                RUNNING_pixels = pixels;
                RUNNING_ang = ang;
            end

        end
        
        global_color_image = DrawWeightImage(RUNNING_pixels,RUNNING_max_rSq);
        
        global_image = RgbBlend(image, global_color_image, 'green');

        fprintf('Maximum angle is %0.2f.\n',RUNNING_ang);
        
    end
    
    
    function [rSq,pixels,ring_lengths] = AnalyseGlobal(edgelist,EDGE_DATA,max_angle)
        
        rSq = NaN * zeros(array_size);
        pixels(array_size(1),array_size(2),array_size(3)).locations = [];
        ring_lengths = zeros(array_size);

        % Ang always positive, so...
        alt_max_ang = max_angle - pi;
        max_angles_dir = [max_angle alt_max_ang];

        % For each edge...

        for e = edge_numbers

            if size(edgelist{e},1) <= 10
                continue
            end

            for direction = 1:length(probe_directions)

                dir_rSq = NaN * zeros(2,length(probe_distances));

                for distance = 1:length(probe_distances)

                    points = EDGE_DATA(e).ShiftedPoints(:,distance,direction);
                    intensities = EDGE_DATA(e).Intensity(:,distance,direction);

                    intensities = intensities(points>0);
                    angles = EDGE_DATA(e).Orientation(points>0);
                    points = points(points>0);

                    pixels(e,direction,distance).locations = points;
                    ring_lengths(e,direction,distance) = size(points,1);

                    if length(points) >= 10

                        for ang = 1:2

                            curr_max_ang_dir = max_angles_dir(ang);

                            current_ang_dst = abs(curr_max_ang_dir - angles);
                            pi_check = current_ang_dst > pi;

                            current_ang_dst = ~pi_check .* current_ang_dst + ...
                                pi_check .* (2*pi - current_ang_dst);

                            [O,ind] = sort(current_ang_dst);
                            G = intensities(ind);
                            P = points(ind);
                            
                            % Check matrix rank
                            X = [O ones(size(P))];
                            [~,R,~] = qr(X,0);
                            p = sum(abs(diag(R)) > max(length(P),2)*eps(R(1)));

                            if p >= 2                            
                                [B,~,~,~,stats] = regress(G,X);
                            else
                                B = 0;
                                stats = 0;
                            end

                            % If slope is above threshold...
                            if (stats(1) > Rsquare_threshold && abs(B(1)) > slope_threshold)
                                dir_rSq(ang,distance) = stats(1);
                            else
                                dir_rSq(ang,distance) = NaN;
                            end

                        end
                    end
                end

                % Use distance 1 to determine which angle was best.

                [~,max_dst1] = max(dir_rSq(:,1));
                rSq(e,direction,:) = dir_rSq(max_dst1,:);

            end

        end
        
    end
    
    function image = DrawWeightImage(pixels,weights)

        % Take structure of point list [pixels] and color each list
        % using appropriate [weights]. Each weight corresponds to a
        % list.

        if size(pixels) ~= size(weights)
            error('Cannot draw weighted image.');
        end

        image = zeros(imsize);

        for i = 1:size(weights,1)
            for j = 1:size(weights,2)
                for k = 1:size(weights,3)
                    p = pixels(i,j,k).locations;
                    w = weights(i,j,k);
                    if w > Rsquare_threshold
                        color = (w-Rsquare_threshold)/(1.5-Rsquare_threshold);
                    else
                        color = 0;
                    end
                    image(p) = color;
                end
            end
        end
    end
    
  
end
    

function new_image = RgbBlend(old_image,alpha_image,channel_string)
   
    % Takes input grey image and amps up a particular color channel
    % according to an alpha image.
    
    new_image = zeros(size(old_image,1),size(old_image,2),3);
    
    if strcmp(channel_string,'red')
        channel = 1;
    elseif strcmp(channel_string,'green')
        channel = 2;
    elseif strcmp(channel_string,'blue')
        channel = 3;
    end
    
    for c = 1:3
        new_image(:,:,c) = (1-alpha_image) .* old_image;
    end
    
    new_image(:,:,channel) = new_image(:,:,channel) + (alpha_image).*(2^16-1);
   
    new_image = uint16(new_image);
    
end



function avg = wmean(x,w,dim)
   
    if nargin < 1
        dim = 1;
    end
    
    % If dim is 0, sum whole array using weights, ignoring NaNs.
    
    if dim == 0
        w = w(~isnan(x));
        x = x(~isnan(x));
        k = w ./ sum(w);
        avg = sum(x.*k);
    else
        k = w ./ sum(w,dim);
        avg = sum(x.*k,dim);
    end
    
end


function color_vector = ang2rgb(angle,value)
    
    % Convert angle in radians [0 pi] to RGB hue with value.
    
    if angle < 0
        angle = angle + pi;
    end
    
    hue = angle / pi;
    color_vector = hsv2rgb([hue 1 value]);
    
end


function [mask_edges,mask_vectors] = AnalyseMaskEdges(mask_binary,mask_vectors)
   
    if (nargin<2 || isempty(mask_vectors))
        
        mask_edges = edge(mask_binary);
        mask180 = skeletonOrientation(mask_edges,[9 9]);
        mask_vectors = deg2rad(GetMaskOrientations(mask_binary,mask180));
        
    end
    
    % Put maskVectors inside [-pi pi]
    mask_vectors(mask_vectors>pi) = mask_vectors(mask_vectors>pi) - 2*pi;
    
end
