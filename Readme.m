% NOTES FOR THE MODEL

% This model aims to classify images comprised only of smooth intensity
% gradients and sharp contours.

% 1. Blur detector. 

% Low-passes the image, and checks if the result is extremely close to the 
% original image. This indicates a lack of high spatial frequencies. 
% (Alternatively, analyse the Fourier frequency component.)

% 2. Self-occlusion detector. 

% High-passes the image to extract sharp edges. The high-pass threshold is
% not extremely high, because even blurry (soft) edges can appear as
% defocused self-occlusions.

% Checks local patches along the edge for linear falloff from the intensity
% maximum as a function of edge orientation. Independent analyses for both
% sides of the edge. Luminance maxima are preferred to NOT be associated
% with "under-edges" (i.e. orientations of about -90).

% 2a. Illumination detector.

% If a self-occlusion is detected, the illumination azimuth is set to the
% luminance maximum along the edge. Only some of the patches will happen to
% capture a luminance maximum - these patches will detect linear falloff on
% BOTH sides of those maximum.

% 2b. Second blur detector.

% A second pass of the blur detector, discounting sharp contours that from
% (2) do not appear as self-occlusions. These contours are locally blurred
% to eliminate them from consideration.

% 3. Curvature detector.

% Takes the self-occlusions detected in (2) and estimates how rapidly the
% surface curves toward frontoparallel as a function of distance from the
% edge.

% Performs the same analysis as (2) at increasing distance D from the
% self-occlusion, and checks to see how the slope of the linear intensity
% falloff from the intensity maxima changes as a function of distance. The
% rate of change of this slope (of slopes) is a cue for changing surface 
% orientation.

% If the falloff slope decreases, then surface slant is decreasing
% (positive curvature - spherical).

% If the falloff slope increases, then surface slant is increasing
% (negative curvature). In practice, this means that the slope next to the
% edge is small, which generally makes the edge more difficult to detect,
% and the surface will also never (smoothly) approach frontoparallel.

% If the falloff slope is constant, then surface slant is not changing
% (conical).