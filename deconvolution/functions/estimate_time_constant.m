function g = estimate_time_constant(y, p, sn, lags, fudge_factor)
%% Estimate noise standard deviation and AR coefficients if they are not present

%% inputs:
%   y: N X T matrix, fluorescence trace
%   p: positive integer, order of AR system
%   sn: scalar, noise standard deviation, estimated if not provided
%   lags: positive integer, number of additional lags where he autocovariance is computed
%    fudge_factor: float (0< fudge_factor <= 1) shrinkage factor to reduce bias

%% outputs
%   g: 1 x p vector, AR coefficient

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% adapted from the MATLAB implemention by Eftychios Pnevmatikakis and the
% Python implementation from Johannes Friedrich

%% References
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments
if ~exist('p', 'var') || isempty(p)
    p = 2;
end
if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(y);
end
if ~exist('lags', 'var') || isempty(lags)
    lags = 5;
end
if ~exist('fudge_factor', 'var') || isempty(fudge_factor)
    fudge_factor = 1;
end

%% estimate time constants 
lags = lags + p;
if ~isempty(which('xcov')) %signal processing toolbox
    xc = xcov(y,lags,'biased');
else
    ynormed = (y - mean(y));
    xc = nan(lags + 1, 1);
    for k = 0:lags
        xc(k + 1) = ynormed(1 + k:end)' * ynormed(1:end - k);
    end
    xc = [flipud(xc(2:end)); xc] / numel(y);
end
xc = xc(:);
A = toeplitz(xc(lags+(1:lags)),xc(lags+(1:p))) - sn^2*eye(lags,p);
g = pinv(A)*xc(lags+2:end);

% while max(abs(roots([1,-g(:)']))>1) && p < 3
%     warning('No stable AR(%i) model found. Checking for AR(%i) model \n',p,p+1);
%     p = p + 1;
%     g = estimate_time_constants(y,p,sn,lags);
% end
% if p == 5
%     g = 0;
% end

% re-adjust time constant values
rg = roots([1;-g(:)]);
if ~isreal(rg); rg = real(rg) + .001*randn(size(rg)); end
rg(rg>1) = 0.95 + 0.001*randn(size(rg(rg>1)));
rg(rg<0) = 0.15 + 0.001*randn(size(rg(rg<0)));
pg = poly(fudge_factor*rg);
g = -pg(2:end);

