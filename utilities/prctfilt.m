function Y = prctfilt(X,p,window,shift)

% percentile filtering along the last dimension
% INPUTS
% X:        input data
% p:        percentile
% window:   window over which to compute the percentile
% shift:    length of window shifting (default: window)

% OUTPUT
% Y:        filtered version

sizX = size(X);
if ndims(X) == 1
    X = X(:)';
elseif ~ismatrix(X)
    X = reshape(X,[],sizX(end));
end

Xpad = [X(:,1:ceil((window-1)/2)),X,X(:,sizX(end) - (ceil((window-1)/2)-1:-1:0))];

blocks_end = window:shift:size(Xpad,2);
blocks_start =  1:shift:size(Xpad,2);
ln = length(blocks_end);
blocks_start(ln+1:end) = [];

Xpcs = zeros(size(Xpad,1),ln);
for i = 1:ln
    Xpcs(:,i) = prctile(Xpad(:,blocks_start(i):blocks_end(i)),p,2);
end

dx = diff(Xpcs,1,2);

Xf = kron(Xpcs(:,1:end-1),ones(1,shift)) + kron(dx,1/shift*(0:shift-1));
Xf = padarray(Xf,[0,size(X,2)-shift*(ln-1)],'post','symmetric');
Y = X - Xf;
Y = reshape(Y,sizX);