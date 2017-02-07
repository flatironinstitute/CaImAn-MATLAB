function [Yest, results] = lle(Y, ssub, rr, ACTIVE_PX)
%% approximate the background with locally-linear embedding
%% inputs:
%   Y: d1*d2*T 3D matrix, video data
%   rr: scalar, average neuron size
%   ssub: spatial downsampling factor
%   ACTIVE_PX:  indicators of pixels to be approximated

%% outputs:
%   Yest: d1*d2*T 3D matrix, reconstructed video data
%   results: struct variable {weights, ssub}
%       weights: d1*d2 cell, each element is a 2*J matrix. Row 1 has the indice of the
%       ring neighbors and row 2 has the corresponding weights. 
%       ssub:    scalar, spatial downsampling factor 

%% Author: Pengcheng Zhou, Carnegie Mellon University,2016

%% input arguments
[d1, d2, T] = size(Y);

% center the fluorescence intensity by its mean
Ymean = mean(Y, 3);
Y = Y - bsxfun(@minus, Ymean, ones(1, 1, T));

% average neuron size
if ~exist('rr', 'var')|| isempty(rr)
    rr = 15;
end
% spatial downsampling
if ~exist('ssub', 'var') || isempty(ssub)
    ssub = 1;
end

%downsample the data
if ssub>1
    Y = imresize(Y, 1./ssub);
    [d1s, d2s, ~] = size(Y);
    rr = round(rr/ssub)+1;
else
    d1s = d1;
    d2s = d2;
end

% pixels to be approximated
if exist('ACTIVE_PX', 'var') && ~isempty(ACTIVE_PX)
    ACTIVE_PX = reshape(double(ACTIVE_PX), d1, d2);
    ACTIVE_PX = (imresize(ACTIVE_PX, 1/ssub)>0);
    ind_px = find(ACTIVE_PX(:));  % pixels to be approxiamted
else
    ind_px = (1:(d1s*d2s));
end

%% determine neibours of each pixel
c_shift = (-rr:rr);
r_shift = round(sqrt(rr^2-c_shift.^2));
c_shift = [c_shift, c_shift(2:end)];
r_shift = [-r_shift, r_shift(2:end)];

[csub, rsub] = meshgrid(1:d2s, 1:d1s);
csub = reshape(csub, [], 1);
rsub = reshape(rsub, [], 1);
csub = bsxfun(@plus, csub, c_shift);
rsub = bsxfun(@plus, rsub, r_shift);
% remove neighbors that are out of boundary
ind = or(or(csub<1, csub>d2s), or(rsub<1, rsub>d1s));
csub(ind) = nan;
rsub(ind) = nan;
% options = optimoptions('lsqlin','Algorithm','active-set', 'Display', 'none');

%% run approximation
 gamma = 0.001; % add regularization
Y = reshape(Y, d1s*d2s, []);
Yest = zeros(size(Y));
weights = cell(d1s, d2s); 
for m=1:length(ind_px)
    px = ind_px(m);
    ind_nhood = sub2ind([d1s,d2s], rsub(px, :), csub(px, :));
    ind_nhood(isnan(ind_nhood)) = [];
    J = length(ind_nhood);
    
    G = bsxfun(@minus, Y(ind_nhood, 1:3:end), Y(px, 1:3:end));
    w = (G*G'+gamma*sum(G(:).^2)*eye(J))\ones(length(ind_nhood), 1);
    w = w/sum(w);
%     w = lsqlin(Y(ind_nhood,:)', Y(px,:)', [], [], ones(1, J), ...
%         1, ones(J, 1)*0, ones(J, 1)*3/J, [], options); % weights have box constraints [min_w, max_w]
    Yest(px, :) = w'*Y(ind_nhood, :);
    weights{px} = [ind_nhood; w']; 
end
results.weights = weights; 
results.ssub = ssub; 

ind = 1:(d1s*d2s);
ind(ind_px) = [];
if ~isempty(ind)
    temp = imfilter(Y, ones(3, 3)/9, 'replicate');
    Yest(ind, :) = temp(ind, :); % without approximation
end
Yest = reshape(Yest, d1s, d2s, []);

    %% return the result
if ssub>1 %up sampling
    Yest = imresize(Yest, [d1, d2]);
end

clear Y;
Ybaseline = Ymean - median(Yest, 3);
Yest = bsxfun(@plus, Yest, Ybaseline);
