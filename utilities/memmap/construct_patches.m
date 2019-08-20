function patches = construct_patches(sizY,siz,overlap,min_size)

% constructs the coordinates of the different patches for batch processing
% INPUTS
% sizY:     size of input data (excluding time dimension)        
% siz:      size of each patch  (default: [32,32,4], if third spatial dimension is present)
% overlap:  amount of overlap in each direction (default: [4,4,2], if third spatial dimension is present)
% min_size: minimum patch size is each direction, otherwise merge with previous patch (defalt: [4,4,2])

% OUTPUT
% patches:  cell structure with start/end points for each patch

% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2016

dimY = length(sizY);    % dimension of dataset (2d or 3d)

if nargin < 4 || isempty(min_size)
    min_size = [8,8,3*ones(1,dimY-2)];
end

if nargin < 3 || isempty(overlap)
    overlap = [4,4,2*ones(1,dimY-2)];
end

if nargin < 2 || isempty(siz)
    siz = [32,32,4*ones(1,dimY-2)];
end

if any(siz<=overlap); error('Size of patch must be greater than the amount of overlap'); end

x_start = 1:siz(1)-overlap(1):sizY(1)-overlap(1);
x_end   = min(x_start + siz(1) - 1,sizY(1));
if (x_end(end) - x_start(end) + 1 < min_size(1)) || (length(x_end) > 1 && (x_end(end) - x_end(end-1) < min_size(1)))    
    x_start(end) = [];
    x_end(end-1) = [];
end

y_start = 1:siz(2)-overlap(2):sizY(2)-overlap(2);
y_end   = min(y_start + siz(2) - 1,sizY(2));
if (y_end(end) - y_start(end) + 1 < min_size(2)) || (length(y_end) > 1 && (y_end(end) - y_end(end-1) < min_size(2)))    
    y_start(end) = [];
    y_end(end-1) = [];
end

if dimY == 3
    z_start = 1:siz(3)-overlap(3):sizY(3)-overlap(3);
    z_end   = min(z_start + siz(3) - 1,sizY(3));
    if z_end(end) - z_start(end) + 1 < min_size(3)
        z_start(end) = [];
        z_end(end-1) = [];
    end    
else
    z_start = 1;
    z_end   = 1;
end

[X1,Y1,Z1] = meshgrid(x_start,y_start,z_start);
[X2,Y2,Z2] = meshgrid(x_end,y_end,z_end);

if dimY == 2
    Z1 = [];
    Z2 = [];
end

patches = mat2cell([X1(:),X2(:),Y1(:),Y2(:),Z1(:),Z2(:)],ones(numel(X1),1),2*dimY);