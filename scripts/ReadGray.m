function grayImage = ReadGray(imagename,color)
    
    if nargin < 2
        color = 0;
    end

% Reads imagename and converts to grayscale.

image = imread(imagename);

switch size(image,3)
    case 1
        grayImage = image;
    case 2
        error('NxNx2 image not understood.');
    case 3
        if color
            grayImage = image;
        else
            grayImage = rgb2gray(image);
        end
    case 4 
        if color
            grayImage = image(:,:,1:3);
        else
            grayImage = rgb2gray(image(:,:,1:3));
        end
end
    
