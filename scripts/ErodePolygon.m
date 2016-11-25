function [eroded_verts,angles,vert_im,eroded_im] = ErodePolygon(verts,erosion,res)
    
    if nargin < 3
        res = 512;
    end
    
    % Fit polygon into 1024/1024 binary image.
    
    min_x = min(verts(:,1));
    min_y = min(verts(:,2));
    max_x = max(verts(:,1));
    max_y = max(verts(:,2));
    
    % Rescale to integer 3-2046
    
    if max_x - min_x >= max_y - min_y
        scale = max_x-min_x;
    else
        scale = max_y-min_y;
    end
    
    s_verts = [verts(:,1)-min_x, verts(:,2)-min_y];
    
    s_verts = round( s_verts * (res-20) / scale + 10 );
    
    % Convert to polygon
    
    vert_im = poly2mask(s_verts(:,1),s_verts(:,2),res,res);
    
    se = strel('disk',round(erosion));
    
    eroded_im = imerode(vert_im,se);
    
    eroded_verts = bwboundaries(eroded_im);
    eroded_verts = eroded_verts{1};
    eroded_verts = [eroded_verts(:,2) eroded_verts(:,1)]; % Convert RC format to XY
    eroded_verts = ((eroded_verts-10) * scale / (res-20)); % Scale back to original coordinates
    eroded_verts = [eroded_verts(:,1)+min_x, eroded_verts(:,2)+min_y];
    
    angles = zeros(size(eroded_verts,1),1);

    angles(1) = atan2(eroded_verts(end,2)-eroded_verts(1,2),eroded_verts(1,1)-eroded_verts(end,1));

    for v = 2:size(eroded_verts,1);

        angles(v) = atan2(eroded_verts(v-1,2)-eroded_verts(v,2),eroded_verts(v,1)-eroded_verts(v-1,1));

    end