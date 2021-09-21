function [IMF,F] = fpFIF2(img,opts)
% Function that decomposes the input image on several IMFs:
% img = IMF(:,:,1) + IMF(:,:,2) + ... + IMF(:, :,end)
% where the last component is the trend.
%
% Function is dedicated to work with optical fringe pattern images. It
% requires 'prefixed_double_filter.mat' presence in working directory.
%
% Inputs:
%   img - input image
%   opts - decomposition options (OPTIONAL VARIABLE)
%   We recommend to leave default options, which should work fine
%   independent of the given fringe pattern. However, for specialized
%   applications opts may be adjusted:
%       opts.alpha - percentile of the detected signal frequencies that
%           are stored in current IMF. Range: 1-100, default 30. Lower the
%           alpha value to obtain more dense frequency sampling
%       opts.Xi - parameter for adjusting the distances between extrema to
%           the mask length. Range - around 1-3, default 1.8. If the
%           image decomposition is unsatisfactory a small Xi change may
%           help
%       opts.p - cosine fraction of Tukey window for additional IMFs
%           filtering. Range - 0-0.9, default 0.6. The larger is p, the
%           smaller amount of low frequencies is stored in additional IMF.
%           Set opts.disp = 1 to see which IMFs are additional
%       opts.maxIMFs - maximal number of IMFs that may be generated
%           (default 50)
%       opts.maxInner - maximal number of iterations in the inner loop
%           (default 200)
%       opts.delta - tolerance parameter to stop the inner loop
%           (default 0.001)
%       opts.minExtr - number of detected extrema that will stop the
%           algorithm (default 3)
%       opts.disp - choose if want to display IMF information in command 
%           window (0 - No, 1 - Yes), default 0
%   Outputs
%       IMF - produced IMFs
%       F - filters used to filter each IMF
%
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland
%
% Last modified: 16.09.2021
%
% Ref. Mikołaj Rogalski, Mateusz Pielach, Antonio Cicone, Piotr Zdańkowski,
%      Luiza Stanaszek, Katarzyna Drela, Krzysztof Patorski, Barbara
%      Łukomska, Maciej Trusiak. "Tailoring 2D Fast Iterative Filtering
%      algorithm for low-contrast optical fringe pattern preprocessing."
%      2021. Submitted
%
% Algorithm is based on FIF2 code: https://github.com/Acicone/FIF2

%% Deal with the input
if nargin == 0
    error('There should be at least one input variable');
elseif nargin == 1
    opts = [];
elseif nargin > 2
    error('There should be max 2 input variables');
end

if isfield(opts,'alpha') == 0; opts.alpha = 30; end
if isfield(opts,'Xi') == 0; opts.Xi = 1.8; end
if isfield(opts,'p') == 0; opts.p = 0.6; end
if isfield(opts,'maxIMFs') == 0; opts.maxIMFs = 50; end
if isfield(opts,'maxInner') == 0; opts.maxInner = 200; end
if isfield(opts,'delta') == 0; opts.delta = 0.001; end
if isfield(opts,'minExtr') == 0; opts.minExtr = 3; end
if isfield(opts,'disp') == 0; opts.disp = 0; end
%% Initialization

% Parameters for third stop criterium
stop3 = sum(img-min(img,[],'all'),'all');
tol3 = 0.05;

S = size(img);
IMF = [];

% Core to generate the filter w_n
load('prefixed_double_filter','MM');

I_n = img;  % Currently processed image

%% fpFIF2 outer loop
n = 1;  % IMF number
m_n = 0;    % Mask length
epsilon_last = 1; % Previous average distance between I_n extrema (pix)
while n < opts.maxIMFs
    
    % Calculating average distance between extrema (epsilon)
    [epsilon, stop12] = CalcEpsilon(I_n, epsilon_last, opts);
    if stop12 == 1; break; end % If stop criterion 1 or 2 is reached
    
    m_nb = m_n; % Previous mask length
    m_n = round(2*opts.Xi*epsilon); % Mask length
    
    if m_n < 1.1*m_nb
        % If new mask length is close to or smaller than previous
        if opts.disp == 1
            disp(['IMF ',num2str(n),' m increased from ', num2str(m_n),...
                ' to ',num2str(ceil(m_n*1.1))]);
        end
        m_n=ceil(m_n*1.1);
        epsilon = m_n/2/opts.Xi;
    end
    
    % Generate the filter
    w_n = get_mask_2D_v3(MM,m_n);
    
    % Precompute the FFT of the signal and of the filter, so we can
    % apply it with almost no computations
    [Sy, Sx] = size(I_n);
    [Swy, Swx] = size(w_n);
    
    % Note that if the filter has a support larger than our image
    % we might need to extend the image a little bit. The repmat
    % command takes care of that.
    if Swy > Sy || Swx > Sx
        ExtL=ceil(max(Swy/Sy, Sx/Swx));
        FI_n = repmat(I_n,ExtL);
        [Sy, Sx] = size(FI_n);
    else
        ExtL = 1;
        FI_n = I_n;
    end
    
    % Pad w_n the right way -- this is required to make sure the
    % center of gravity of the filter is in position (1,1) before
    % taking the FFT
    l1 = floor(Swy/2); m1 = floor(Swx/2);
    l2 = Swy - l1; m2 = Swx - m1;
    
    Fw_n = zeros(Sy, Sx);
    Fw_n(1:l2,1:m2) = w_n(l1+1:end,m1+1:end);
    Fw_n(end-l1+1:end,end-m1+1:end) = w_n(1:l1,1:m1);
    Fw_n(end-l1+1:end,1:m2) = w_n(1:l1,m1+1:end);
    Fw_n(1:l2,end-m1+1:end) = w_n(l1+1:end,1:m1);
    
    % Precomputing FFTs for later use
    FI_n = fft2(FI_n);
    Fw_n = fft2(Fw_n);
    
    % r1 and r2 are updated throughout the iterations so that r1 /
    % r2 is equal to the relative change in the solution at step j.
    % To accomplish these, some values are accumulated in
    % filter_factors (which indeed describes the action of the
    % accumulated filter on all the frequencies), making use of the
    % variable incr_filter.
    r2 = abs(FI_n).^2;
    r1 = r2 .* abs(Fw_n).^2;
    incr_filter = (1 - abs(Fw_n)).^2;
    filter_factors = ones(size(Fw_n));
    
    
    f_n = S(2)/epsilon/2;   % Detected image frequency
    f_n_last = S(2)/epsilon_last/2; % Previous frequency
    
    % Additional IMF generation - high frequency case
    if epsilon < 8
        % Additional filter
        Fw_nT = FreqMask(FI_n,f_n,f_n_last,opts,1);
        
        if ~isempty(Fw_nT)
            FI_n = FI_n .* Fw_nT;   % Updating the FI_n
            upd = real(ifft2(FI_n));
            IMF(:,:,n) = I_n - upd(1:Sy/ExtL,1:Sx/ExtL); % Generating IMF
            I_n = I_n - IMF(:,:,n); % Updating the I_n
            if opts.disp == 1; disp(['IMF ',num2str(n),...
                    ' additional filtering - high frequency case']); end
            F{n} = Fw_nT;
            n = n+1;
            
        end
    end
    % Additional IMF generation - low frequency case
    if epsilon >= 8 && epsilon/epsilon_last >= 1.8 && epsilon < 100
        Fw_nT = FreqMask(FI_n,f_n,f_n_last,opts,2);
        if ~isempty(Fw_nT)
            %             Flt(:,:,n) = Fw_nT;
            FI_n = FI_n .* Fw_nT;   % Updating the FI_n
            upd = real(ifft2(FI_n));
            IMF(:,:,n) = I_n - upd(1:Sy/ExtL,1:Sx/ExtL); % Generating IMF
            I_n = I_n - IMF(:,:,n); % Updating the I_n
            if opts.disp == 1; disp(['IMF ',num2str(n),...
                    ' additional filtering - low frequency case']); end
            F{n} = Fw_nT;
            n = n+1;
        end
    end
    if n>=opts.maxIMFs
        if opts.disp == 1; disp('Stop - max number of IMFs reached'); end
        break
    end
    epsilon_last = epsilon;
    
    %% inner loop
    inStepN = 0;
    SD=1;
    while SD>opts.delta && inStepN < opts.maxInner
        inStepN=inStepN+1;
        
        % Construct the residual for checking the stopping
        % criterion for inner loop
        h_avenrm = sqrt(sum(sum( r2 .* filter_factors )));
        SD = sum(sum( r1 .* filter_factors )) ./ ...
            h_avenrm^2;
        
        if SD <= opts.delta || inStepN >= opts.maxInner
            Fw_n = (1-Fw_n).^(inStepN);
            FI_n = FI_n.*(1-Fw_n);  % Updating the FI_n
            upd = real(ifft2(FI_n));
            IMF(:,:,n) = I_n - upd(1:Sy/ExtL,1:Sx/ExtL); % Generating IMF
            I_n = I_n - IMF(:,:,n); % Updating the I_n
            if opts.disp == 1; disp(['IMF ',num2str(n),' m = ',...
                    num2str(m_n), ' inner step = ', num2str(inStepN)]); end
            F{n} = 1-Fw_n;
            n = n+1;
        else
            % Update the filter factors for the next round
            filter_factors = incr_filter .* filter_factors;
        end
    end
    
    if n>=opts.maxIMFs
        if opts.disp == 1; disp('Stop - max number of IMFs reached'); end
        break
    end
    
    % Third stop criterion
    if sum(I_n - min(I_n,[],'all'),'all') < tol3*stop3
        if opts.disp == 1
            disp('Stop - third stop criterion (background consist almost no signal)');
        end
        break
    end
end
% last IMF is the residual
IMF(:,:,end+1) = I_n;
% fl = 1-Fw_n;
end

%% Auxiliary functions
function A = get_mask_2D_v3(w,k)
% get the mask with length 2*k+1 x 2*k+1
% k must be integer
% w is the area under the curve for each bar
% A  the mask with length 2*k+1 x 2*k+1

L=length(w);
m=(L-1)/2;  %2*m+1 =L punti nel filtro
w=[w zeros(1,(L-1)/2)];
A=zeros(k+1,k+1);
if k<=m % The prefixed filter contains enough points
    if mod(k,1)==0     % if the mask_length is an integer
        for i=0:k
            for j=0:k
                d1=sqrt(i^2+j^2);
                d=k;
                if d1 > d
                    A(j+1,i+1)=0;
                else
                    s=(m-1)+L/2*d1/d;
                    t=s+2;
                    s2=ceil(s)-s;
                    t1=t-floor(t);
                    A(j+1,i+1)=sum(w(ceil(s):floor(t)))+s2*w(ceil(s))+t1*w(floor(t));
                end
            end
        end
        A=[rot90(A(2:end,2:end),2) rot90(A(2:end,:),3)';rot90(A(:,2:end),1)' A];
        A=A/sum(sum(A,1),2);
    else   % if the mask length is not an integer
        disp('Need to write the code!')
        A=[];
        return
    end
else % We need a filter with more points than MM, we use interpolation
    disp('Need to write the code!')
    A=[];
    return
end

end

function Fw_nT = FreqMask(FI_n,f_n,f_n_last,opts,mode)
% Function that creates Fourier spectrum amplitude filter to generate the
% additional IMFs
%   Inputs:
%       FI_n - image Fourier spectrum
%       f_n - image frequency
%       f_n_last - previous image frequency
%       opts - decomposition options
%       mode - filter option (1 - high frequency case, 2 - low frequency
%           case)
%   Output:
%       Fw_nT - generated filter

% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland

tmp0 = 0;
[Sy, Sx] = size(FI_n);
[xx,yy] = meshgrid(-round(Sx/2):round(Sx/2)-1,-round(Sx/2):Sx/Sy:round(Sx/2)-1);
ff = sqrt(xx.^2+yy.^2); % Frequencies in image

if mode == 1
    % mode 1 - high frequency case
    tol = 0.05; % Tolarance parameter - to allow generating the Fw_nT
    mask = zeros(Sy, Sx);
    mask(ff>f_n) = 1; mask(ff>f_n_last) = 0;
    % Image frequencies between f_n_last and f_n
    frqces = mask.*abs(fftshift(FI_n));
    [fy,fx] = find(frqces == max(frqces(:)));
    f_0 = ff(fy(1),fx(1));  % Found maximal frequency
    
    fqDiff = f_n_last - f_n;
    if f_0 - f_n > tol*fqDiff && f_n_last - f_0 > tol*fqDiff
        tmp0 = 1;
    end
    
elseif mode == 2
    % mode 2 - low frequency case
    f_0 = f_n;
    tmp0 = 1;
end

% Generating the filter basing on f_0 frequency
if tmp0 == 1
    tmp = tukeywin(1000,opts.p);
    fq = sum(tmp==1)/1000/2;
    
    % 1D Tukey window - 1D filter template
    T1D = tukeywin(round(f_0/fq),opts.p);
    qq = length(T1D);
    if qq<Sx
        extra = zeros(round((Sx-qq)/2),1);
        T1D = [extra;T1D;extra];
    end
    T1D = T1D(ceil(end/2):end);
    
    % Generating 2D Tukey window
    [x,y] = meshgrid(-Sx/2:Sx/2-1,-Sx/2:Sx/Sy:Sx/2-1);
    r = round(sqrt(x.^2+y.^2))+1;
    Fw_nT = zeros(Sy,Sx);
    for tt = 1:Sy
        for ss = 1:Sx
            hh = r(tt,ss);
            if hh <= Sx/2
                Fw_nT(tt,ss) = T1D(hh);
            end
        end
    end
    
    % Generated filter
    Fw_nT = ifftshift(Fw_nT);
else
    Fw_nT = [];
end

end

function [epsilon, stop12] = CalcEpsilon(I_n, epsilon_last, opts)
% Function that returns the average distance between extrema (epsilon)
%   Inputs:
%       I_n - image
%       epsilon_last - epsilon calculated in previous iteration (to ensure
%           that epsilon will be larger than epsilon_last)
%       opts - decomposition options
%   Outputs:
%       epsilon - average distance between extrema
%       stop12 - equals 1 if first or second fpFIF2 stop criterion is
%           fulfilled. If not then stop12 = 0
%
% Created by:
%   Mikołaj Rogalski,
%   mikolaj.rogalski.dokt@pw.edu.pl
%   Institute of Micromechanics and Photonics,
%   Warsaw University of Technology, 02-525 Warsaw, Poland

stop12 = 0;

[nY,nX] = size(I_n);

% I_n preprocessing
in = I_n - imgaussfilt(I_n,1);

% Calculating the extrema map (E_n)
msk = strel('disk',1);
MaxDil = imdilate(in,msk); %MAX filtration
MinDil = -imdilate(-in,msk); %MIN filtration
MaxMap = ~(in - MaxDil); %binary map of maxima
MinMap = ~(in - MinDil); %binary map of minima
E_n = MaxMap + MinMap;
E_n(1:end,1) = 0;E_n(1:end,end) = 0;
E_n(1,1:end) = 0;E_n(end,1:end) = 0;

% Number of extremes
NoE = sum(sum(E_n));

if NoE <= opts.minExtr
    % If first stop criterion is fulfilled
    stop12 = 1;
    epsilon = [];
    if opts.disp == 1
        disp('Stop - first stop criterion (number of extrema smaller than or equal to opts.minExtr)');
    end
else
    % Distances map
    D_n = bwdist(E_n);
    % figure; imagesc(D)
    % Finding the maxima (M_n) of the D_n
    msk = strel('square',3);
    MaxDil = imdilate(D_n,msk); %MAX filtration
    M_n = ~(D_n - MaxDil); %binary map of maxima
    
    % Map of distances between extrema
    DBE_n = 2*D_n.*M_n;
    DBE_n(DBE_n<epsilon_last) = 0;
    
    % Dividing the DBE_n on 100 regions and calculating local epsilon for 
    % each region
    si = round(min(nY,nX)/10);
    L_epsilon = zeros(1,floor(nY/si)*floor(nX/si));
    uu = 0;
    for m2 = 1:si:nX-si+1
        for m1 = 1:si:nY-si+1
            uu = uu+1;
            L_epsilon(uu) = mode(double(nonzeros(DBE_n(m1:m1+si-1,m2:m2+si-1))));
        end
    end
    
    % Final epsilon value
    epsilon = prctile(L_epsilon,opts.alpha);
    
    if round(2*opts.Xi*epsilon) > max(nY,nX)
        % If second stop criterion is fulfilled
        stop12 = 1;
        if opts.disp == 1
            disp('Stop - second stop criterion (mask size larger than the image size)');
        end
    end
end
end