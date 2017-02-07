function sn = GetSn(Y, range_ff, method)
%% Estimate noise standard deviation

%% inputs:
%   Y: N X T matrix, fluorescence trace
%   range_ff : 1 x 2 vector, nonnegative, max value <= 0.5, range of frequency (x Nyquist rate) over which the spectrum is averaged
%   method: string, method of averaging: Mean, median, exponentiated mean of logvalues (default)

%% outputs:
%   sn: scalar, std of the noise

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% adapted from the MATLAB implemention by Eftychios Pnevmatikakis and the
% Python implementation from Johannes Friedrich

%% References
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments
if ~exist('range_ff', 'var') || isempty(range_ff)
    range_ff = [.25, .5];
end
if ~exist('method', 'var') || isempty(method)
    method = 'logmexp';
end
if any(size(Y)==1)
    Y = reshape(Y, [], 1);
else
    Y = Y';
end

%% estimate the noise
[psdx, ff] = pwelch(Y, [],[],[], 1);
indf = and(ff>=range_ff(1), ff<=range_ff(2));
switch method
    case 'mean'
        sn=sqrt(mean(psdx(indf, :)/2));
    case 'median'
        sn=sqrt(median(psdx(indf,:)/2));
    case 'logmexp'
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));    
    otherwise
        fprintf('wrong method! use logmexp instead.\n'); 
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));
end
sn = sn';




