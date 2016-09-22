function Cn = correlation_image_3D(Y,sz,sizY)

% Construct correlation image based on neighboring pixels. Works for both
% 2D and 3D imaging data
% 
% INPUTS:
% Y: raw data either in (N-D format or as a 2-D matrix)
% sz: define the relative location of neighbours. it can be scalar, 2
%       element vector or a binary matrix
%       scalar: number of nearest neighbours
%                   For 2D datasets either 4 or 8, default 4 
%                   For 3D datasets either 6 or 26, default 6
%       2 element vector: [rmin, rmax], the range of neighbours. the
%           distance between the neighbour and the pixel is d, dmin <= r < dmax.
%       matrix: a binary matrix indicating the location of neighbors
% sizY: spatial dimensions of data
%
% OUTPUT:
% Cn: correlation image

% Authors: Eftychios A. Pnevmatikakis, Simons Foundation, 2016
% and Pengcheng Zhou, Carnegie Mellon University, 2016


% reshape the data
if nargin<3 && ismatrix(Y)
    error('Data has to be in N-d format or dimensions need to be provided');
elseif ismatrix(Y)
    Y = reshape(Y,[sizY,numel(Y)/prod(sizY)]);
else
    sizY = size(Y);
    sizY(end) = [];
end

dimY = length(size(Y))-1;
sizY = sizY(1:dimY);

if (nargin<2) || isempty(sz)
    if dimY == 2
        sz = 4;
    else
        sz = 6;
    end
end

% construct neighboring filter
if islogical(sz)
    assert(ndims(sz)==dimY,'Binary matrix must have the same number of spatial dimensions as the data.')
    SZ = double(sz);
elseif isscalar(sz)    
    if sz == 4         
        SZ = [0,1,0; 1,0,1; 0,1,0]; 
    elseif sz == 8
        SZ = [1,1,1; 1,0,1; 1,1,1];
    elseif sz == 6
        SZ(:,:,1) = [0,0,0; 0,1,0; 0,0,0]; 
        SZ(:,:,2) = [0,1,0; 1,0,1; 0,1,0];
        SZ(:,:,3) = [0,0,0; 0,1,0; 0,0,0];
    elseif sz == 26
        SZ = ones(3,3,3);
        SZ(2,2,2) = 0;
    end
else
    % the specified neighbours has a distance within the domain [dmin,
    % dmax)
    sz = ceil(sz);
    dmin = min(sz); dmax = max(sz);
    xsub = (-dmax+1):(dmax-1);      % x subscript
    ysub = rsub;                    % y subscript
    if dimY == 2 
        [xind, yind] = meshgrid(xsub, ysub);
        R = sqrt(xind.^2+yind.^2);
        SZ = (R>=dmin) .* (R<dmax);
    else
        [xind, yind, zind] = meshgrid(xsub, ysub, zsub);
        R = sqrt(xind.^2+yind.^2+zind.^2);
        SZ = (R>=dmin) .* (R<dmax);
    end
end

% normalize the data
mY = mean(Y,dimY+1);
Y = bsxfun(@minus,Y,mY);
sY = sqrt(mean(Y.*Y,dimY+1));
Y = bsxfun(@times,Y,1./sY);

% compute correlation image
Yconv = imfilter(Y, SZ);        % sum over the neighbouring pixels
MASK = imfilter(ones(sizY), SZ);   % count the number of neighbouring pixels
Cn = mean(Yconv.*Y, dimY+1)./MASK;   % compute correlation and normalize