function [Orientations] = skeletonOrientation(skel,blksz,angthreshold,edgesubs)
% SKELETONORIENTATION Calculate the local orientation of a skeleton
%
%Inputs: 
%  skel:  MxN binary skeleton image.
%
%  blksz: Size of block to look around for local orientation
%         Both elements must be odd and greater than or equal to three
%         blksz can be a 1x2 vector meaning [row x col] block size or 1x1
%               for square block 
%         optional, defaults to [5 5]
%
%Outputs:
%  Orientation: image that is the size of skel with zeros outside of skel
%               and orientations elsewhere
%
%

%
% Copyright 2013 The MathWorks, Inc. 
% SCd 9/25/2013
%

    %Error checking:
    assert(nargin>=1||nargin<=3,'One to three inputs expected');
    assert(islogical(skel),'skel should be logical');
    assert(ismatrix(skel),'skel should be a matrix');    
    sz = size(skel);
    assert(isequal(sz,size(skel)),'Sizes of M and C expected to be equal');
    
    %Handling blksz
    
    if nargin <= 2
        angthreshold = 0;
    end
        
    if nargin == 1
        %Default
        blksz = [5 5];
    else
        %assertions
        assert(all(blksz>=3),'blksz elements must be greater than or equal to three');
        assert(all(mod(blksz,2)==1),'blksz elements expected to be odd');
        if isscalar(blksz)
            blksz = blksz([1,1]);
        else
            assert(isequal(size(blksz),[1 2]),'blksz is expected to be a variable of size [1x1] or [1x2]');
        end
    end
    
    %Find the skeleton pixels' index
    if nargin < 4
        [row,col] = find(skel);
    else
        row = edgesubs(:,1);
        col = edgesubs(:,2);
    end
    npts      = numel(row);
    
    %Pad the array and offset the rows/cols so every local block fits
    padAmount = floor(blksz./2); %distance from center to edge of block
    skelPad   = padarray(skel,padAmount); %We need to pad so that image boundary pixels are contained in a block
    
    %Preallocate Orientations
    Orientations = zeros(sz); 
    
    %Some parameters
    %-Bottom of block will be the same as center before pad
    %-Top will be bottom + block size - 1 (inclusive    
    rowHigh = row+blksz(1)-1; 
    colHigh = col+blksz(2)-1;
    center  = padAmount+1; %Center of small block
    
    %Now the engine, we will loop over the image, create each local block
    %of pixels that are touching the index of interest. Do a connected
    %components analysis, remove pixels not connected to the center pixel
    %and then calculate orientation.  
    
    %Start the engine!
    for ii = 1:npts;
        %Extract small block
        block = skelPad(row(ii):rowHigh(ii),col(ii):colHigh(ii));
        
        %Label and calculate orientation
        Label = bwlabel(block);
        center_label = Label==Label(center(1),center(2)); %only label of center pixel
        rp    = regionprops(center_label,'Orientation');
        
        %Set orientation of the center pixel equal to the calculated one
        
        if (ii > 1 && angthreshold)
        
            % Check if orientation or (orientation + pi) is closer to
            % previous angle and choose closer one
            
            prevAng = Orientations(row(ii-1),col(ii-1));
            angA = rp.Orientation;
            angB = rp.Orientation + 180;
            
            % Bring angB into [-180 180) range.
            if angB >= 180
                angB = angB - 360;
            elseif angB < -180
                angB = angB + 360;
            end
            %angB = angB
            
            % Compare distances to prevAng
            dstA = abs(prevAng - angA);
            if dstA > 180
                dstA = 360 - dstA;
            end
            dstB = abs(prevAng - angB);
            if dstB > 180
                dstB = 360 - dstB;
            end
            
            % Set new orientation to closest value
            if dstA < dstB
                %display('Choosing angA');
                Orientations(row(ii),col(ii)) = angA;
            elseif dstB < dstA
                %display('Choosing angB');
                Orientations(row(ii),col(ii)) = angB;                
            end

            % DEBUG
%             figure(1)
%             imshow(Orientations,[-180 180]);
%             WaitSecs(0.1);

        else
            Orientations(row(ii),col(ii)) = rp.Orientation;    
        end
    end
    
end