function [g, sn] = estimate_parameters(fluor, p, range_ff, method, lags, fudge_factor)
%% Estimate noise standard deviation and AR coefficients if they are not present 

%% inputs: 
%   fluor: N X T matrix, fluorescence trace 
%   p: positive integer, order of AR system
%   lags: positive integer, number of additional lags where he autocovariance is computed
%   range_ff : 1 x 2 vector, nonnegative, max value <= 0.5, range of frequency (x Nyquist rate) over which the spectrum is averaged
%   method: string, method of averaging: Mean, median, exponentiated mean of logvalues (default)
%    fudge_factor: float (0< fudge_factor <= 1) shrinkage factor to reduce bias

%% outputs 
%   g: 1 x p vector, AR coefficient 
%   sn: scalar, std of the noise 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments 
if ~exist('p', 'var') || isempty(p)
    p = 2; 
end 
if ~exist('range_ff', 'var') || isempty(range_ff)
    range_ff = [.25, .5]; 
end 
if ~exist('method', 'var') || isempty(method)
    method = 'logmexp'; 
end 
if ~exist('lags', 'var') || isempty(lags)
    lags = 5; 
end 
if ~exist('fudge_factor', 'var') || isempty(fudge_factor)
    fudge_factor = 1.0; 
end 

%% estimate noise level 
fluor = reshape(fluor, 1, []); 
sn = GetSn(fluor, range_ff, method); 

%% estmate AR coefficients 
g = estimate_time_constant(fluor, p, sn, lags, fudge_factor); 