function finalImage = SpatialNoise2(imdim,fpower,fSD,numberSamples,aniso,gammaP,seed)

record = 1;
    
if nargin < 4
    numberSamples = 100;
end

if nargin < 6
    gammaP = 1;
end

if nargin < 7
    rng('shuffle');
else
    rng(seed);
end

if nargin < 5
    aniso = 0;
end

if nargin < 2
    fpower = 10*rand;
end

if nargin < 3
    fSD = 0.25;
end

% Choose band(s) of permitted frequency powers. Min frequency is 1 (power
% zero). Max frequency is 1024 (power 10). PDF is Gaussian over frequency
% powers.

freqPowerMean = fpower
freqPowerSD = fSD

% Add up numberSamples of sine waves with random frequencies, orientations, and phases

finalImage = zeros(imdim,imdim);

[Xm, Ym] = meshgrid(1:imdim,1:imdim);

% Define orientation density function

x = 1:629;
if aniso ~= 0
    sd = 50;
    centre1 = 0.25*629;
    centre2 = 0.75*629;
    mode1 = gauss(x,1,centre1,sd/aniso);
    mode2 = gauss(x,1,centre2,sd/aniso);
    th_pdf = mode1 + mode2;
    th_pdf = th_pdf - min(th_pdf); % Kill horizontal orientations.
else
    th_pdf = x;
end
%th_pdf = 1/(2*pi) + (butterfly/(2*pi)) * cos((x-629/4)/50);

x = 0:0.01:4;
th_samples = randp(th_pdf,1,numberSamples);
figure(8), plot(x,gauss(x,1,freqPowerMean,freqPowerSD)), axis([0 4 0 1]);
title('Frequency distribution');
figure(10), plot(th_pdf);
title('Orientation distribution');
axis([1 629 0 1.1]);
%histogram(th_samples);

frequencies = zeros(numberSamples,1);
orientations = zeros(numberSamples,1);

for i = 1:numberSamples
    
    freq = 10 ^ normrnd(freqPowerMean,freqPowerSD) / imdim;        % Random frequency
    frequencies(i) = freq;
    amp = 1/freq;
    phi = rand*2*pi;                                              % Random phase
    theta = th_samples(i)/100;                                    % Random orientation
    orientations(i) = theta;
    
    Xt = Xm * cos(theta);
    Yt = Ym * sin(theta);
    XYt = Xt + Yt;
    XYf = XYt * freq * 2*pi;
    grating = amp*sin(XYf + phi);
    
    finalImage = finalImage + grating;
    
end

finalImage = finalImage / numberSamples;

% Gamma adjustment
finalImage = screenGamma(finalImage,gammaP);

minI = min(finalImage(:));
maxI = max(finalImage(:));

range = [minI-(maxI-minI)/7, maxI];

%figure(9), scatter(orientations,frequencies);
figure(11), imshow(finalImage,range);
axis off; axis image;

%WriteHDR(finalImage,strcat('HDR_sine_image_',num2str(fpower),'_',num2str(fSD)));
%WriteImage(finalImage,strcat('sine_image_',num2str(fpower),'_',num2str(fSD)),'tif');
if record
    Write16BitImage(finalImage,strcat('sine_image_',num2str(fpower),'_',num2str(fSD)),'tif');
end

end

function image_out = screenGamma(image_in,power)
    
    x = 0:0.001:1;
    y = x.^power;
    %figure(12), plot(x,y);
    image_out = image_in.^power;
    
end

function f = gauss(x,a,b,c)
    
    % Gaussian with height a, mean b, and standard deviation c
    
    exponent = -(x-b).^2 / (2*c^2);
    f = a*exp(exponent);
    
end