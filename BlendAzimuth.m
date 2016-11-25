function blended_image = BlendAzimuth(image_name,freq,save_name)
    
    % Blends multiple azimuth renders together according to random noise.
    
    if nargin < 3
        save_name = 'temp';
    end
    if nargin < 2
        freq = 0.5;
    end
    
    res = 1024;
    
    blend_noise = SpatialNoise2(res,freq,0.1,1000,0,1);
    blend_noise = uint8(90*(blend_noise - min(blend_noise(:)))/(max(blend_noise(:))-min(blend_noise(:))));

    figure(1);
    imshow(blend_noise,[0 90]);
    
    % Black means use azimuth 0, white means use azimuth 180
    
    images = zeros(res,res,91);
    
    for im = 0:4:360
        
        number = num2str(im);
        
        if im < 100
            number = strcat('0',number);
        end
        if im < 10
            number = strcat('0',number);
        end
        
        index = im/4 + 1;
        image = ReadGray(sprintf('images/azimuth/%s/%s_azimuth0%s.tif',image_name,image_name,number));
        images(:,:,index) = image;
        
    end
    
    blended_image = zeros(res,res);
    
    for row = 1:res
        for col = 1:res
            blended_image(row,col) = images(row,col,blend_noise(row,col)+1);
        end
    end
    
    blended_image = uint16(blended_image);
    figure(2);
    imshow(blended_image);
    
    try
        imwrite(blended_image,sprintf('blended_image%s.tif',save_name),'tif');
    catch
        display('Couldn''t write image.');
    end
end