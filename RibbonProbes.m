function RibbonProbes
    
    % Create random interpolated ribbon image within parameters specified 
    % below, then run gauge figure task along all probes on a particular
    % strip of that ribbon.
    %
    % Only relevant ribbon parameters are accessed here. Vertices, image 
    % size, etc. are all set in the RibbonShader code.
    
    %% Parameters
    
    % Gauge figure parameters
    
    number_probes = 16;             % Number of probes per image
    probe_radius = 5;               % Size of probe disc and rod
    probe_resolution = 10;          % Resolution of probe
    probe_color = [255 0 0];        % Color of probe
    
    
    % Ribbon parameters
    
    image_size =    256;
    
    offset = 0;   %   [x1 x2 x3 x4 ; ...
                  %   y1 y2 y3 y4];
          
    amplitude = 0;%   [x1 x2 x3 x4 ; ...
                  %   y1 y2 y3 y4];
                 
    azimuth = 0;  %   [x1 x2 x3 x4 ; ...
                  %   y1 y2 y3 y4];
                 
    background =    .5;
    
    falloff =       1;
    
    
    % IDW parameters
    
    IDWweight = [-0.5 -1 -1.5 -2 -2.5 -3];
    IDWnoise = 0:0.25:1;
    
    %% Create IDW
    
    % Create a random contour for every combination of IDWweight and
    % IDWnoise. Currently 20 images.
    
    set = 1;
    num = 0;
    
    for w = IDWweight
        for n = IDWnoise
            num = num + 1;
            DetectSelfOcclusions(set,w,n,num2str(num));
            %CreateIDW(set,w,n,num2str(num));
            while ~KbCheck;
                KbCheck;
            end
        end
    end
    
end


function DetectSelfOcclusions(set,IDWweight,IDWnoise,num)
    
    image_name = sprintf('set%d/IDW_weight%0.2f_noise%0.2f_contour%s',set,IDWweight,IDWnoise,num);
    image = ReadGray(strcat('images/ribbon_images/',image_name,'.tif'));
    
    SelfOcclusionDetector(image,0.05);
    WaitSecs(1);
    
end
    

function [image, probes] = CreateIDW(set,IDWweight,IDWnoise,num)
    
    image_name = sprintf('set%d/IDW_weight%0.2f_noise%0.2f_contour%s',set,IDWweight,IDWnoise,num);

    RibbonShader('Set',set,'Name',num,'IDW_Weight',IDWweight,'Noise',IDWnoise);

    image = ReadGray(strcat('images/ribbon_images/',image_name,'.tif'));
    load(strcat('images/ribbon_images/',image_name,'_probes.mat'));

    probes = edgelink(edge(double(probe_binary)));

end