function Y = prctfilt(X,p,window,shift,mode)

% percentile filtering along the last dimension
% INPUTS
% X:        input data
% p:        percentile
% window:   window over which to compute the percentile
% shift:    length of window shifting (default: window)
% mode:     if mode == 1 then returns Y - baseline, otherwise returns baseline

% OUTPUT
% Y:        filtered version

if nargin < 2 || isempty(p); p = 20; end
if nargin < 3 || isempty(window); window = 200; end
if nargin < 4 || isempty(shift); shift = 200; end
if nargin < 5 || isempty(mode); mode = 1; end

sizX = size(X);
if ndims(X) == 1
    X = X(:)';
elseif ~ismatrix(X)
    X = reshape(X,[],sizX(end));
end

window = min(window,sizX(end));
shift = min(shift,window);

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
if mode == 1
    Y = X - Xf;
else
    Y = Xf;
end
Y = reshape(Y,sizX);