function contrasted_image = ContrastDrop(image,mask,varargin)
   
    if isa(image,'char')
        image = ReadGray(sprintf('images/%s.tif',image));
    end 
    if isa(mask,'char')
        mask = ReadGray(sprintf('images/masks/%s.tif',mask));
    end
    
    if isa(image,'uint16')
        bit_depth = 16;
    elseif isa (image,'uint8')
        bit_depth = 8;
    else
        disp('Image is not an integer array.');
    end
    
    gamma = 1;
    maxdist = [];
    precrop = [1 1 size(image,1) size(image,2)];
    postcrop = precrop;
    global_contrast = 1;
    contrast_minmax = [0 1];
    CIRCLE = 0;
    
    for arg = 1:length(varargin)
        if strcmp(varargin{arg},'Gamma')
            gamma = varargin{arg+1};
        elseif strcmp(varargin{arg},'Precrop')
            precrop = varargin{arg+1};
        elseif strcmp(varargin{arg},'Postcrop')
            postcrop = varargin{arg+1};
        elseif strcmp(varargin{arg},'Global')
            global_contrast = varargin{arg+1};
        elseif strcmp(varargin{arg},'CircleUp')
            CIRCLE = 1;
            circle_pos = varargin{arg+1};
            circle_radius = varargin{arg+2};
            circle_minmax = varargin{arg+3};
            circle_gamma = varargin{arg+4};
        elseif strcmp(varargin{arg},'CircleDown')
            CIRCLE = -1;
            circle_pos = varargin{arg+1};
            circle_radius = varargin{arg+2};
            circle_minmax = varargin{arg+3};
            circle_gamma = varargin{arg+4};
        end
    end
    
    % PRE-CROP
    image = image(precrop(2):precrop(4),precrop(1):precrop(3),:);
    mask = mask(precrop(2):precrop(4),precrop(1):precrop(3),:);
    
    mask = double(mask);
    mask = mask / max(mask(:));
    
    image_d = double(image);
    image_d = 2^bit_depth * (image_d - min(image_d(:))) / (max(image_d(:)) - min(image_d(:)));
    
    binary = imbinarize(mask);
    dists = double(bwdist(binary));
    if isempty(maxdist)
        maxdist = max(dists(:));
    end
    dists = min(1, dists / maxdist);
    %dists = (1 - dists);
    dists = dists .^ gamma;
    dists = dists * (contrast_minmax(2)-contrast_minmax(1)) + contrast_minmax(1);
    dists = imgaussfilt(dists,2);
    
    %dists = min(dists,2);
    
    % Mean luminance
    meanlum = mean(image_d(:));
    
    BG = meanlum;
    %BG = 0;
    
    % Transform luminance mean to zero
    centered_image = image_d - meanlum;
    
    % As dist grows, contrast goes down (gamma function)
    contrasted_image = (centered_image .* dists .* global_contrast);
    
    % OPTIONAL: Circular contrast manipulation
    if CIRCLE~=0
        circle_filter = zeros(size(contrasted_image));
        circle_filter(circle_pos(1),circle_pos(2)) = 1;
        circle_filter = bwdist(circle_filter);
        circle_filter = min(1, circle_filter/circle_radius) .^ circle_gamma;
        if CIRCLE==1
            circle_filter = 1 - circle_filter;
        end        
        circle_filter = circle_filter * (circle_minmax(2)-circle_minmax(1)) + circle_minmax(1);
        contrasted_image = contrasted_image .* circle_filter;
        
    end
    
    % Re-center mean
    contrasted_image = contrasted_image + meanlum;
    
    % OPTIONAL: Scale down contrasted image if it is outside range
%     contrasted_image = 2^bit_depth * (contrasted_image - min(contrasted_image(:))) ...
%         / (max(contrasted_image(:)) - min(contrasted_image(:)));
    
    % POST-CROP
%     image_d = image_d(postcrop(2):postcrop(4),postcrop(1):postcrop(3),:);
%     contrasted_image = contrasted_image(postcrop(2):postcrop(4),postcrop(1):postcrop(3),:);
%     mask = mask(postcrop(2):postcrop(4),postcrop(1):postcrop(3),:);

    % Mask images
    image_d = (mask)*meanlum + (1-mask).*image_d;
    contrasted_image = (mask)*BG + (1-mask).*contrasted_image;
    
    if bit_depth == 8
        image_d = uint8(image_d);
        contrasted_image = uint8(contrasted_image);
    elseif bit_depth == 16
        image_d = uint16(image_d);
        contrasted_image = uint16(contrasted_image);
    end
    
    figure(1);
    imshow(image_d);
    figure(2);
    imshow(contrasted_image);
    
end