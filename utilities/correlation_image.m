function Cn = correlation_image(Y,sz,d1,d2,flag_norm)

% construct correlation image based on neighboing pixels
% Y: raw data
% sz: define the relative location of neighbours. it can be scalar, 2
%       element vector or a binary matrix
%       scalar: number of nearest neighbours, either 4 or 8, default 4
%       2 element vector: [rmin, rmax], the range of neighbours. the
%           distance between the neighbour and the pixel is d, dmin <= r <
%           dmax.
%       matrix: a squared matrix (2k+1)*(2k+1)indicating the location of neighbours
% d1,d2: spatial dimensions
% flag_norm: indicate whether Y has been normalized and centered ( 1 is
%   yes, 0 is no)

% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2015
% with modifications from Pengcheng Zhou, Carnegie Mellon University, 2015

%% preprocess the raw data
if nargin<5;     flag_norm = false; end
if ismatrix(Y)
    Y = reshape(Y,d1,d2,size(Y,2));
end
if ~strcmpi(class(Y),'single'); Y = single(Y); end
if (nargin<2);   sz = [0,1,0; 1,0,1; 0,1,0]; end

[d1,d2,~] = size(Y);    % image dimension

if (~flag_norm)
    % centering
    mY = mean(Y,3);
    Y = bsxfun(@minus, Y, mY);
    % normalizing
    sY = sqrt(mean(Y.*Y, 3));
    Y = bsxfun(@times, Y, 1./sY);
end

%% construct a matrix indicating location of the matrix
if  isscalar(sz)
    if sz == 8      % 8 nearest neighbours
        sz = [1,1,1; 1,0,1; 1,1,1];
    elseif sz==4
        sz = [0,1,0; 1,0,1; 0,1,0];
    end
elseif length(sz(:)) == 2
    % the specified neighbours has a distance within the domain [dmin,
    % dmax)
    sz = ceil(sz);
    dmin = min(sz); dmax = max(sz);
    rsub = (-dmax+1):(dmax-1);      % row subscript
    csub = rsub;      % column subscript
    [cind, rind] = meshgrid(csub, rsub);
    R = sqrt(cind.^2+rind.^2);
    sz = (R>=dmin) .* (R<dmax);
end

%% compute the correlation
Yconv = imfilter(Y, sz);        % sum over the neighbouring pixels
MASK = imfilter(ones(d1,d2), sz);   % count the number of neighbouring pixels
Cn = mean(Yconv.*Y, 3)./MASK;   % compute correlation and normalize
