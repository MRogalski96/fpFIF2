%% Input fringe pattern image

imS = 1000; % Image size
sig1x = 3;   % Gaussian background sigma parameter (x direction)
sig1y = 3;   % Gaussian background sigma parameter (y direction)
sig2 = 0.3; % Noise variance
D = 2;      % Phase dynamic range
Tx = 80;    % Carrier fringes x period
Ty = 60;    % Carrier fringes y period

backgroundGT = 5 + gausswin(imS,sig1y)*gausswin(imS,sig1x)';
noiseGT = sig2*randn(imS);
[x,y] = meshgrid(1:imS);
phsGT = peaks(imS)*D;   % Ground truth phase
fringesGT = cos(phsGT + (x/Tx+y/Ty)*2*pi);    % Ground truth fringe pattern
image = backgroundGT+fringesGT+noiseGT;   % input image 

figure;imagesc(image); colormap gray; axis image; title('Input image')

%% fpFIF2 

% Default version
[IMFs,F] = fpFIF2(image);

% % With custom options:
% options = [];
% options.disp = 1;
% options.alpha = 20;
% IMF = fpFIF2(image,options);

%% Displaying IMFs
figure;
cn = ceil(sqrt(size(IMFs,3)+1));
cn2 = ceil((size(IMFs,3)+1)/cn);
subplot(cn2,cn,1); imagesc(image); colormap gray; colorbar; title('Input image'); %axis image
for tt = 2:size(IMFs,3)+1
    subplot(cn2,cn,tt); imagesc(IMFs(:,:,tt-1)); colormap gray; colorbar;
    title(['fpFIF2; IMF ',num2str(tt-1),'/',num2str(size(IMFs,3))]); %axis image
end

%% Results
modes = 6:14;   % IMFs that consist fringes
fringes = sum(IMFs(:,:,modes),3);
noise = sum(IMFs(:,:,1:(modes(1)-1)),3);
background = sum(IMFs(:,:,(modes(end)+1):end),3);

figure;
subplot(2,4,1); imagesc(image); axis image; title('Input image'); colorbar
subplot(2,4,2); imagesc(noiseGT); axis image; 
title('Ground truth noise'); colorbar
subplot(2,4,3); imagesc(fringesGT,[-1 1]); axis image; 
title('Ground truth fringes'); colorbar
subplot(2,4,4); imagesc(backgroundGT); axis image; 
title('Ground truth background'); colorbar
subplot(2,4,6); imagesc(noise,[min(noiseGT(:)) max(noiseGT(:))]); 
axis image; title('Filtered noise'); colorbar
subplot(2,4,7); imagesc(fringes,[-1 1]); axis image; 
title('Filtered fringes'); colorbar
subplot(2,4,8); imagesc(background,[min(backgroundGT(:)) max(backgroundGT(:))]); 
axis image; title('Filtered background'); colorbar
colormap gray

%% fpFIF2 for timelapse
% use F returned for one image to decompose second, simiar image
fringesGT2 = cos(phsGT + (-x/Ty+y/Tx)*2*pi); % grounth truth fringes
image2 = backgroundGT+fringesGT2+noiseGT;   % input image 

[noise2,fringes2,background2] = fpFIF2_for_timelapse(image2,F,[modes(1),modes(end)]);

%% Results

figure;
subplot(2,4,1); imagesc(image2); axis image; title('Input image'); colorbar
subplot(2,4,2); imagesc(noiseGT); axis image; 
title('Ground truth noise'); colorbar
subplot(2,4,3); imagesc(fringesGT2,[-1 1]); axis image; 
title('Ground truth fringes'); colorbar
subplot(2,4,4); imagesc(backgroundGT); axis image; 
title('Ground truth background'); colorbar
subplot(2,4,6); imagesc(noise2,[min(noiseGT(:)) max(noiseGT(:))]); 
axis image; title('Filtered noise'); colorbar
subplot(2,4,7); imagesc(fringes2,[-1 1]); axis image; 
title('Filtered fringes'); colorbar
subplot(2,4,8); imagesc(background2,[min(backgroundGT(:)) max(backgroundGT(:))]); 
axis image; title('Filtered background'); colorbar
colormap gray