function BlurDetector(image)
   
    % Calculates high-pass image and detects edges. If not enough edges can
    % be detected (sum of edge array), then the image is considered blurry.
    
    % If enough edges are detected, we move to phase two: checking if they
    % are self-occlusions that could cause the low-frequency image
    % components to appear focused. (SelfOcclusionDetector)
    
    % Any edges without significant correlations are erased (with some
    % padding) from the high pass image (i.e. set to 0) and edges are
    % searched for again.
    
    % Currently, this method will report that a flat grey square will be
    % "blurry." In order to detect if there is some luminance variation
    % (devoid of 
    
end 


function ImageFourier(image)

    F = fft2(image);
    F = fftshift(F);

    F = abs(F);
    F = log(F+1);
    F = mat2gray(F);

    figure(3);
    imshow(F,[]);

end