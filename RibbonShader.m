function [INTERP_IMAGE, probe_binary] = RibbonShader(varargin)
    
    % Create a random ribbon contour and shade it according to specified
    % amplitude, offset, and azimuth parameters.
    
    NOISE = 0;
    shading_mode = 'Bezier';
    lum_method = 'Falloff';
    name = '';
    set = '';
    
    % Options
    
    SAVE_VECTORS = 0;                   % Don't save vectors if random.
    DO_IDW = 1;                         % Output IDW interpolation.
    RECORD = 0;                         % Record image to file.
    ANALYSE = 0;                        % Analyse gradients using mask.

    
    % Bezier contour parameters
    
    res = 512;
    n_vertices = 8;
    radius_minmax = [0.5 1.5];
    
    
    % Ribbon parameters
    
    ribbon_width = 0;
    n_ribbons = 1;
    probe_distance = 20; % Which of the n_ribbons to set as probe origins?
    ribbon_gap = round(ribbon_width / n_ribbons);
    
    
    % Bezier ribbon parameters
    
    offset(1,:) = [0 .33 .66 1] * (n_ribbons-1);
    offset(2,:) = [.5 .5 .5 .5];
    
    amplitude(1,:) = [0 .33 .66 1] * (n_ribbons-1);
    amplitude(2,:) = [.5 .5 0 0];
    
    l_azimuth(1,:) = [0 .33 .66 1] * (n_ribbons-1);
    l_azimuth(2,:) = -pi/2 * ones(1,4);
    
    l_elevation = pi/4;
    l_intensity = 1;
    
    background_color = .5;
    
    falloff = 1;
    
    
    % IDW parameters
    
    if ribbon_gap > 0        
        sample_radius = 1.3*ribbon_gap;
    else
        sample_radius = res;
    end
    d_weight = -1;
    
    
    % Display options
    
    fade_ribbons = 0;
    gamma = 1;
    crop = 0;
    range = [0 1];

    
    % Try loading pre-existing ribbon vectors
    try
        load(strcat('data/RibbonVectors/',num2str(r),'_res',num2str(res),'.mat'))
    catch
        display('Could not find ANY ribbon vectors. Creating all...')
        Binary = zeros(res,res,n_ribbons);
        Edges = zeros(res,res,n_ribbons);
        Distances = zeros(res,res,n_ribbons);
        Vectors = zeros(res,res,n_ribbons);
    end
    

    for arg = 1:length(varargin)        
        if strcmp(varargin{arg},'Azimuth')
            l_azimuth = varargin{arg+1};
        elseif strcmp(varargin{arg},'Name')
            name = varargin{arg+1};
        elseif strcmp(varargin{arg},'Set')
            set = varargin{arg+1};
        elseif strcmp(varargin{arg},'Noise')
            NOISE = varargin{arg+1};
        elseif strcmp(varargin{arg},'Offset')
            offset = varargin{arg+1};
        elseif strcmp(varargin{arg},'Amplitude')
            amplitude = varargin{arg+1};
        elseif strcmp(varargin{arg},'Elevation')
            l_elevation = varargin{arg+1};
        elseif strcmp(varargin{arg},'Intensity')
            l_intensity = varargin{arg+1};
        elseif strcmp(varargin{arg},'IDW_Weight')
            d_weight = varargin{arg+1};
        elseif strcmp(varargin{arg},'IDW_Radius')
            sample_radius = varargin{arg+1};
        elseif strcmp(varargin{arg},'Res')
            res = varargin{arg+1};
        elseif strcmp(varargin{arg},'Vertices')
            n_vertices = varargin{arg+1};
        elseif strcmp(varargin{arg},'Radius')
            radius_minmax = varargin{arg+1};
        elseif strcmp(varargin{arg},'Record')
            RECORD = 1;
        elseif strcmp(varargin{arg},'Analyse')
            ANALYSE = 1;
        elseif strcmp(varargin{arg},'Gamma')
            gamma = varargin{arg+1};
        elseif strcmp(varargin{arg},'Crop')
            crop = varargin{arg+1};
        elseif strcmp(varargin{arg},'Range')
            range = varargin{arg+1};
        end

    end

    % Generate outer contour    
    verts = RandomBoundary(res,n_vertices,radius_minmax);
    
    % Specify amplitude and offset weights
    
    if strcmp(shading_mode,'Lambertian')
        
        sim_slant_angle = ((n_ribbons-1):-1:0) * (pi/2) / (n_ribbons-1);
        offset_points = cos(sim_slant_angle) * cos(l_elevation);
        amplitude_points = sin(sim_slant_angle) * sin(l_elevation);
        background_color = cos(l_elevation);
        
    elseif strcmp(shading_mode,'Bezier')
        
        offset_points = createBezier(n_ribbons,offset);        
        amplitude_points = createBezier(n_ribbons,amplitude);                
        l_azimuth_points = createBezier(n_ribbons,l_azimuth);

    end
    
    %OFFSET = sin(l_elevation) * sin(s_slant) * ones(1,n_ribbons);
    %AMPLITUDE= cos(l_elevation) * cos(s_slant) * ones(1,n_ribbons);
    
    % Generate N ribbons and create combined image
    EDGE_ROWS = [];
    EDGE_COLS = [];
    EDGE_IMAGE = zeros(res);
    LUM_IMAGE = EDGE_IMAGE;
    
    MASK_INFO = [];
    SAVE_NEW = 0;
    
    for i = 1:n_ribbons
        
        ribbon = (i-1)*ribbon_gap + 1;
        
        if (~isempty(Binary) && size(Binary,3) >= ribbon)
            checkRibbon = Binary(:,:,ribbon);
        else
            checkRibbon = zeros(res);
        end
        
        if ~isequal(checkRibbon,zeros(res))
            currentBinary = Binary(:,:,ribbon);
            currentEdges = Edges(:,:,ribbon);
            currentVectors = Vectors(:,:,ribbon);
            display(strcat(num2str(ribbon),'-th ribbon detected.'));
        else
            SAVE_NEW = 1;
            display(strcat(num2str(ribbon),'-th ribbon could not be found. Creating...'));
            [~,~,outerBinary,currentBinary] = ErodePolygon(verts,ribbon-1,res);
            Binary(:,:,ribbon) = currentBinary;
            currentEdges = edge(currentBinary);
            Edges(:,:,ribbon) = currentEdges;
            currentDistances = bwdist(currentEdges);
            Distances(:,:,ribbon) = currentDistances;
            currentVectors = deg2rad(GetMaskOrientations(currentBinary,skeletonOrientation(currentEdges,[9 9])));
            Vectors(:,:,ribbon) = currentVectors;                        
        end
        
        % If this is the desired probe ribbon, output that binary image. It
        % will later be ordered and used to create spaced probe points.
        if i == probe_distance
            probe_binary = Binary(:,:,ribbon);
        end
        
        if ribbon == 1
            [~,~,outerBinary] = ErodePolygon(verts,1,res);
        end
        
        if strcmp(lum_method,'Cosine')
            currentLum = l_intensity * max(0, offset_points(i) + amplitude_points(i).*cos(currentVectors-l_azimuth_points(i)));
        elseif strcmp(lum_method,'Falloff')
            ang_dst = abs((currentVectors-pi)+l_azimuth_points(i));
            check = ang_dst > pi;
            ang_dst = (~check).*ang_dst + check.*(2*pi - ang_dst);            
            currentLum = l_intensity * max(0, ...
                lerp(offset_points(i)+amplitude_points(i),offset_points(i)-amplitude_points(i),ang_dst/(falloff*pi)));
        end
        
        [currentEdgeRows,currentEdgeCols] = find(currentEdges); % RC coordinates of all edges.
        EDGE_ROWS = [EDGE_ROWS; currentEdgeRows];
        EDGE_COLS = [EDGE_COLS; currentEdgeCols];

        EDGE_IMAGE = ~currentEdges.*EDGE_IMAGE + currentEdges;
        LUM_IMAGE = ~currentEdges.*LUM_IMAGE + currentEdges.*currentLum;

        display(strcat(num2str(i),' ribbons done'));
    end
    
    %% Create probe ring
    
    [~,~,~,probe_binary] = ErodePolygon(verts,probe_distance,res);
        
    if RECORD
           
        bin_name = sprintf('IDW_weight%0.2f_noise%0.2f_contour%s_probes',d_weight,NOISE,name);        
        save(strcat('images/ribbon_images/',bin_name,'.mat'),'probe_binary');

    end 
    
    MASK_INFO.Binary = Binary;
    MASK_INFO.Edges = Edges;
    MASK_INFO.Distances = Distances;
    MASK_INFO.Vectors = Vectors;
    
    if (SAVE_NEW && SAVE_VECTORS)
        save(strcat('data/RibbonVectors/',num2str(r),'_res',num2str(res),'.mat'),'-struct','MASK_INFO');
    end
    
    RIBBON_IMAGE = ~EDGE_IMAGE.*background_color + EDGE_IMAGE.*LUM_IMAGE;
    
    figure(1)
    imshow(RIBBON_IMAGE .^ gamma,[0 1]);
    
    if DO_IDW
        % Interpolate with IDW

        % Get all non-station points within outer boundary:
        if ribbon_width == 0
            [interpRows,interpCols] = find(~EDGE_IMAGE.*outerBinary); % RC coordinates of all strip points.
            figure(10);
            imshow(~EDGE_IMAGE.*outerBinary);            
        else
            [interpRows,interpCols] = find(~EDGE_IMAGE.*outerBinary.*~currentBinary); % RC coordinates of all strip points.
        end 

        lumInd = sub2ind([res res],EDGE_ROWS,EDGE_COLS);
        lumVec = LUM_IMAGE(lumInd);

        INTERP_IMAGE = IDW(EDGE_COLS,EDGE_ROWS,lumVec,interpCols,interpRows,d_weight,'fr',sample_radius,res,1);
        INTERP_IMAGE = INTERP_IMAGE + EDGE_IMAGE.*LUM_IMAGE; % Restore station points

        % Get background mask
        BG_MASK = ~outerBinary .* ~im2bw(INTERP_IMAGE);
        
        % Make anc scale noise
        if NOISE > 0
            noise_image = SpatialNoise2(res,0.3,0.1,1000,0,1);
            noise_image = linscale(noise_image,0.3,1);
        
            INTERP_IMAGE = NOISE.*noise_image + (1-NOISE).*INTERP_IMAGE;
        end
        
        % Add background
        INTERP_IMAGE = BG_MASK * background_color + ~BG_MASK .* INTERP_IMAGE;

        % Create and add alpha for middle
        if fade_ribbons > 0
            ALPHA_MASK = Distances(:,:,ribbon-(fade_ribbons*ribbon_gap))/(fade_ribbons*ribbon_gap) .* Binary(:,:,ribbon-(fade_ribbons*ribbon_gap));
            ALPHA_MASK(ALPHA_MASK>1) = 1;        
        else
            ALPHA_MASK = Binary(:,:,ribbon-(fade_ribbons*ribbon_gap));
        end

        %INTERP_IMAGE = ALPHA_MASK * background_color + (1-ALPHA_MASK) .* INTERP_IMAGE;

        figure(2);
        imshow(INTERP_IMAGE,range);
                
    end
    
    if crop > 0    
        CROPPED_IM = maskCropped.*DISPLAY_IM + (1-maskCropped).*BG;
        figure(3)
        imshow(CROPPED_IM,range);
    end
      
    
    if RECORD
        
        WRITE_IM = uint8(255*INTERP_IMAGE);
        
        im_name = sprintf('set%d/IDW_weight%0.2f_noise%0.2f_contour%s',set,d_weight,NOISE,name);
        
%         im_name = strcat('az',num2str(l_azimuth_range(1)),'-',num2str(l_azimuth_gamma),'-',num2str(l_azimuth_range(2)), ...
%             '_amp',num2str(amp_range(1)),'-',num2str(amp_gamma),'-',num2str(amp_range(2)), ...
%             '_off',num2str(off_range(1)),'-',num2str(off_gamma),'-',num2str(off_range(2)));
        imwrite(WRITE_IM,strcat('images/ribbon_images/',im_name,'.tif'),'tif');
        
    end
    
    if ANALYSE
        
        % For every pixel in the mask...
        LUM_ARRAY = INTERP_IMAGE(maskBinary);
        TILT_ARRAY = tiltArray(maskBinary);
        DIST_ARRAY = maskDistances(maskBinary);
        
        [LUM_LIST,i] = sort(LUM_ARRAY);
        TILT_LIST = TILT_ARRAY(i);
        DIST_LIST = DIST_ARRAY(i);
        
        figure(3)
        scatter3(TILT_LIST,DIST_LIST,LUM_LIST,0.1);
        
        % Now do steps of distances (20)
        for dist = 5:5:80
            DIST_IN = logical((DIST_ARRAY<dist) .* (DIST_ARRAY>=(dist-5)));
            TILT_IN = TILT_ARRAY(DIST_IN);
            LUM_IN = LUM_ARRAY(DIST_IN);
            figure(4);
            subplot(4,4,dist/5);
            scatter(TILT_IN,LUM_IN,0.1);
            axis([-4 4 0 1]);
        end
        
    end
    
    INTERP_IMAGE = uint8(255*INTERP_IMAGE);
    commandwindow;
   
end

function curve = createBezier(n,anchors)
    
    % Fill in value curve using control points for cubic Bezier curve.
    % Takes four input (x,y) anchors and calculates n points on the
    % resulting Bezier curve.
    
    if n <= 1
        curve = anchors(2,:);
        return
    end
    
    x_anchors = anchors(1,:);
    y_anchors = anchors(2,:);

    axis([0 n-1 -2 2]); grid on;
    xlabel u; ylabel Pu;
    
    if length(y_anchors) == 1
        y_anchors = y_anchors * ones(1,4);
    end
    
    % Calculate geometric matrix B:
    
    B = [x_anchors' y_anchors' zeros(4,1)];
    
    % Specify Bezier basis transformation matrix M:
 
    M = [-1  3 -3  1 ; ...
          3 -6  3  0 ; ...
         -3  3  0  0 ; ...
          1  0  0  0 ];
    
    % Calculate coefficient matrix A:
    
    A = M*B;
 
    % Define u-axis
 
    u = linspace(0,1,n);
    u = u';
    unit = ones(size(u));
    U = [u.^3 u.^2 u unit];
    
    % Calculate values for each value of u
 
    Pu = U*A;
    curve = Pu(:,2);
    
    % Plot control polygon (blue) and Bezier curve (red)
 
    line(B(:,1), B(:,2), B(:,3))
    hold on;
    plot3(Pu(:,1), Pu(:,2), Pu(:,3),'r','linewidth',1.0) 
    legend('Polygon','Bezier Curve')
    
end

function point = lerp(a,b,p)
    
    % Interpolate point p in [0 1] between a and b.
    
    point = a + (b-a) * p;
    
end

function array = linscale(array,minV,maxV)
    
    % Scales array into range [minV, maxV]

    arraymin = min(array(:));
    arraymax = max(array(:));

    array = array - arraymin; % Scale to [0 arraymax-arraymin]
    array = array * (maxV - minV) / (arraymax-arraymin); % Scale to [0 maxV-minV]
    array = array + minV; % Scale to [minV maxV]
   

end