function [Cn, PNR] = correlation_image_endoscope(Y, options)
%% compute correlation image of endoscopic data. it has to spatially filter the data first
%% Input:
%   Y:  d X T matrx, imaging data
%   options: struct data of paramters/options
%       d1:     number of rows
%       d2:     number of columns
%       gSiz:   maximum size of a neuron
%       nb:     number of background
%       min_corr: minimum threshold of correlation for segementing neurons
%% Output:
%       Cn:  d1*d2, correlation image
%       PNR: d1*d2, peak to noise ratio
%% Author: Pengcheng Zhou, Carnegie Mellon University. zhoupc1988@gmail.com

%% use correlation to initialize NMF
%% parameters
d1 = options.d1;        % image height
d2 = options.d2;        % image width
gSig = options.gSig;    % width of the gaussian kernel approximating one neuron
gSiz = options.gSiz;    % average size of neurons
sig = 5;    % thresholding noise by sig*std()

if ismatrix(Y); Y = reshape(Y, d1, d2, []); end;  % convert the 3D movie to a matrix
Y(isnan(Y)) = 0;    % remove nan values
Y = double(Y);
T = size(Y, 3);

%% preprocessing data
% create a spatial filter for removing background
psf = fspecial('gaussian', round(gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;

% divide data into multiple patches
patch_sz = [3, 3];
r0_patch = round(linspace(1, d1, 1+patch_sz(1)));
c0_patch = round(linspace(1, d2, 1+patch_sz(2)));
nr_patch = length(r0_patch)-1; 
nc_patch = length(c0_patch)-1; 
Cn = zeros(d1, d2);
PNR = zeros(d1,d2);

% compute correlation_image patch by patch
bd = round(gSiz);
for mr = 1:nr_patch
    r0 = max(1, r0_patch(mr)-bd); % minimum row index of the patch
    r1 = min(d1, r0_patch(mr+1)+bd-1); % maximum row index of the patch
    for mc = 1:nc_patch
        c0 = max(1, c0_patch(mc)-bd); % minimum column index of the patch
        c1 = min(d2, c0_patch(mc+1)+bd-1); % maximum column index of the patch
        
        % take the patch from the raw data
        nrows = (r1-r0+1);  % number of rows in the patch
        ncols = (c1-c0+1);  %number of columns in the patch
        Ypatch = double(Y(r0:r1, c0:c1, :));
        
        % spatially filter the data
        HY = imfilter(Ypatch, psf, 'replicate');
        
        % copute signal to noise ratio
        HY = reshape(HY, [], T);
        HY = bsxfun(@minus, HY, median(HY, 2));
        HY_max = max(HY, [], 2);
        Ysig = get_noise_fft(HY, options);
        tmp_PNR = reshape(HY_max./Ysig, nrows, ncols);
        PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR);
        
        
        %  compute loal correlation
        HY(bsxfun(@lt, HY, Ysig*sig)) = 0;
        tmp_Cn = correlation_image(HY, [1,2], nrows, ncols);
        Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn);
        
    end
end