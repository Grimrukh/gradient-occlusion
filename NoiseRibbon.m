function blended_image = NoiseRibbon(ribbon_image,alpha)
    
    % Assumes 16-bit ribbon image
    
    if isa(ribbon_image,'char')
        load(sprintf('images/%s_pixels.mat',ribbon_image));
        ribbon_image = double(ReadGray(sprintf('images/%s.tif',ribbon_image)))/(2^16-1);
    else
        error('Please give string for ribbon name.');
    end
    
    noise_image = double(ReadGray('sine_image_1.25_0.13.tif'))/(2^16-1);
    noise_image = noise_image(513:1536,513:1536) .* ribbon_pixels;
    
    noise_image = noise_image - 0.5;
    noise_image = 2*noise_image + 0.5;
    
    blended_image = ribbon_image*(1-alpha) + noise_image*(alpha);
    blended_image = blended_image.*ribbon_pixels + ribbon_image.*(1-ribbon_pixels);
    
    % Gamma
    blended_image = blended_image .^ 1.3;
    
    % Put between [0.2 0.8]
    blended_image = blended_image * 0.8 + 0.1;
    
    blended_image = uint16((2^16)*blended_image);
    figure(1);
    imshow(blended_image);
    
end