function [u,v,w] = fpFIF2_for_timelapse(img,F,num)
% Function that performs fast fpFIF2 decomposition of the fringe pattern 
% data for a data series. It requires to firstly perform fpFIF2 on one of
% the images in data series
% 
% Inputs:
%   img - fringe pattern image
%   F - set of filters returned by the fpFIF2 algorithm
%   num - IMFs that consists fringe pattern 
%       num = [n_first, n_last]; 
%           n_first - first IMF that consist fringes
%           n_last - last IMF that consist fringes
% Outputs:
%   u - img noise component
%   v - img fringe component
%   w - img background component
% 
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
% 
% Last modified: 21.09.2021

[Sy,Sx] = size(img);
F1 = 1;
for tt = 1:(num(1)-1)
    if size(F{tt},1) > Sy
        F{tt} = imresize(F{tt},[Sy,Sx]);
    end
    F1 = F{tt}.*F1;
end
spct = fft2(img).*F1;
tmp = real(ifft2(spct));    % fringes + background
u = img-tmp;    % noise
F2 = 1;
for tt = num(1):num(2)
    if size(F{tt},1) > Sy
        F{tt} = imresize(F{tt},[Sy,Sx]);
    end
    F2 = F{tt}.*F2;
end
w = real(ifft2(spct.*F2));  % background
v = tmp - w;    % fringes
end

