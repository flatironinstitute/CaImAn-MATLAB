function Y_ds = dsData(Y, options)
%% downsampling data
% input:
%   Y:      d1*d2*T matrix.
%   time
%   options: struc data containing parameters
%       ssub:   downsampling factor for spatial
%       tsub:   downsampling factor for temporal
% output:
%   Y_ds:   downsampled data with (d1/ssub)*(d2/ssub)*(t/tsub)

% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

% input arguments are not enough
if nargin<2; Y_ds = Y;  return; end

% downsampling factors
if (isfield(options, 'ssub')); ssub = options.ssub; else ssub = 1; end
if (isfield(options, 'tsub')); tsub = options.tsub; else tsub = 1; end
if or(isempty(ssub), ssub<1); ssub = 1; end
if or(isempty(tsub), tsub<1); tsub = 1; end
ssub = round(ssub);
tsub = round(tsub);

% data dimension
[d1, d2, T] = size(Y);

%% downsampe
d1s = ceil(d1/ssub);
d2s = ceil(d2/ssub);
Ts = floor(T/tsub);

if ssub>1
    % spatial downsampling
    Y_ds = imresize(Y, [d1s, d2s], 'box');
else
    Y_ds = Y; 
end

if tsub>1  
    % temporal downsampling, take the mean of continuours tsub frames 
    % if mod(T, tsub)~=0, discard the rest 
    Y_ds = squeeze(mean(reshape(Y_ds(:, :, 1:(Ts*tsub)), d1s, d2s, tsub, Ts), 3)); 
end