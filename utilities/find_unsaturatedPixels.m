function normalPixels = find_unsaturatedPixels(Y, saturationValue, saturationThreshold, saturationTime)

% finds the pixels/voxels that are saturated and returns the ones that are
% not. A pixel is defined as saturated if its observed fluorescence is
% above saturationThreshold*saturationValue at least saturationTime
% fraction of the time.

% Written by Weijian Yang and Eftychios A. Pnevmatikakis, based on an idea
% from Weijian Yang and Darcy Peterka

if nargin < 4 || isempty(saturationTime)
    saturationTime = 0.005;
end

if nargin < 3 || isempty(saturationThreshold)
    saturationThreshold = 0.9;
end

if nargin < 2 || isempty(saturationValue)
    saturationValue = 2^nextpow2(max(Y(:)))-1;
end

Ysat = (Y >= saturationThreshold*saturationValue);
normalPixels = find(mean(Ysat,ndims(Y))<saturationTime);