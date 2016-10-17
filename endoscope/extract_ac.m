function [ai, ci, ind_success] = extract_ac(HY, Y, ind_ctr, sz)
%% given a patch of raw & high-pass filtered calcium imaging data, extract
% spatial and temporal component of one neuron (ai, ci). if succeed, then
% return an indicator ind_succes with value 1; otherwise, 0.
%% inputs:
%       HY:     d X T matrix, filtered patch data
%       Y:      d X T matrix, raw data
%       ind_ctr:        scalar, location of the center
%       sz:         2 X 1 vector, size of the patch

%% Author: Pengcheng Zhou, Carnegie Mellon University.

%% parameters 
nr = sz(1);
nc = sz(2);
min_corr = 0.7;
min_pixels = 5;

%% find pixels highly correlated with the center
HY(HY<0) = 0;       % remove some negative signals from nearby neurons
y0 = HY(ind_ctr, :);
tmp_corr = reshape(corr(y0', HY'), nr, nc);
data = HY(tmp_corr>min_corr, :);


%% estimate ci with the mean or rank-1 NMF
ci = mean(data, 1);
% ci = ci - min(ci); % avoid nonnegative baseline
% [~, ci] = nnmf(ci, 1);
if norm(ci)==0
    ai=[];
    ind_success=false;
    return;
end

%% extract spatial component
% estiamte the background level using the boundary
y_bg = median(Y(tmp_corr(:)<min_corr, :), 1); % using the median of the whole field (except the center area) as background estimation

% sort the data, take the differencing and estiamte ai
thr_noise = 3;      % threshold the nonzero pixels to remove noise
[~, ind_sort] = sort(y_bg, 'ascend');   % order frames to make sure the background levels are close within nearby frames
dY = diff(Y(:, ind_sort), 2, 2);    % take the second order differential to remove the background contributions
dY(bsxfun(@gt, dY, -thr_noise*get_noise_fft(dY))) = 0;   % threshold the result to pick only frames with large signals
dci = diff(ci(ind_sort), 2);
ci_noise = get_noise_fft(dci);
dci(dci>- thr_noise * ci_noise) = 0;
ai = max(0, dY*dci'/(dci*dci'));  % use regression to estimate spatial component

% post-process ai by bwlabel
thr_ratio = 50;    % thresholding nonzero pixels by its ratio with the maximum intensity
temp = full(ai>max(ai)/thr_ratio);
l = bwlabel(reshape(temp, nr, nc), 4);   % remove disconnected components
temp(l~=l(ind_ctr)) = false;
ai(~temp(:)) = 0;
if sum(ai(:)>0) < min_pixels %the ROI is too small
    ind_success=false;
    return;
end

% refine ci again
ind_nonzero = (ai>0);
% ci0 = ci;
ai_mask = mean(ai(ind_nonzero))*ind_nonzero;
ci = (ai-ai_mask)'*ai\((ai-ai_mask)'*Y);
% if get_noise_fft(ci)>get_noise_fft(ci0) % refined ci should have smaller noise
%     ci = ci0;
% end

% set the baseline to be 0
dci = [0, 0, ci(4:end)-ci(1:(end-3))];
ci = ci - median(ci(dci>0));
sn = get_noise_fft(ci);
ci = ci / sn;
ai = ai * sn;
% ci(ci < -thr_noise) = 0;
% return results
if norm(ai)==0
    ind_success= false;
else
    ind_success=true;
end