function maskAngles = GetMaskOrientations(mask,maskOrientations)

% Calculates the orientation of each point on the edges in the input mask
% image. Edges are oriented in 360-degree space. A zero-degree orientation
% is a horizontal edge with black below and white above.

imSize = size(mask);

%% Get orientation up to 180 degrees

maskEdges = edge(mask);

maskAngles = maskOrientations;
maskAngles(maskEdges==0) = NaN; % Set non-edge pixels to NaN

%% Correct orientations by checking mask inside/outside

checkMask = mask;

for row = 1:imSize(1)
    for col = 1:imSize(2)        
        if maskEdges(row,col) == 1
            
            % Check if (x,y) is at the image boundary.
            
            scale = 5;
            
            if (col <= scale && row < imSize(1)-scale) % Left boundary, excluding bottom
                shift = [0 1]; % Check below
            elseif (col >= imSize(2)-scale && row < imSize(1)-scale) % Right boundary, excluding bottom
                shift = [0 1]; % Check below
            elseif (col > scale && col < imSize(2)-scale && row <= scale) % Top boundary, excluding sides
                if maskAngles(row,col) > 0
                    shift = [1 0]; % Check right
                else
                    shift = [-1 0]; % Check left
                end
            elseif (row >= imSize(1)-scale) % Bottom boundary
                if maskAngles(row,col) > 0
                    shift = [1 0]; % Check right
                else
                    shift = [-1 0]; % Check left
                end
            else
                
                if (maskAngles(row,col) <= -45)                    
                    shift = [-1 0];                    
                elseif (maskAngles(row,col) <= 45 && maskAngles(row,col) > -45)                    
                    shift = [0 1];                    
                elseif (maskAngles(row,col) >= 45)
                    shift = [1 0];
                end
                
            end
                        
            shift = shift * scale; % Scale jump up. One pixel is too small.
            
            try
                checkMask(row+shift(2),col+shift(1)) = 0;
            catch
                figure(66);
                imshow(checkMask);
                display(row)
                display(col)
                display(shift)
                error('Bad check.');
            end
            
            if mask(row+shift(2),col+shift(1)) == 1
                maskAngles(row,col) = maskAngles(row,col);
            else
                maskAngles(row,col) = maskAngles(row,col) + 180;
            end
            
        end        
    end
end
