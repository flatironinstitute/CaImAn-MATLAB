function [level,cdf_val] = estimate_percentile_level(data,window,shift)

% Estimates the level of percentile filtering to be used in DF/F extraction
% using a kernel density estimator. 

if ~exist('window','var'); window = 4000; end
if ~exist('shift','var'); shift = 2000; end

data = data(:);
T = length(data);
start_point = 1:shift:T;
end_point = [window:shift:T,T];

min_ln = min(length(start_point),length(end_point));
start_point = start_point(1:min_ln);
end_point = end_point(1:min_ln);
ln_seg = end_point - start_point + 1;

if ln_seg(end) < 2*window/3
    start_point(end) = [];
    end_point(end-1) = [];
    ln_seg = end_point - start_point + 1;
end

level = zeros(length(start_point),1);
cdf_val = zeros(length(start_point),1);

for i = 1:length(start_point)
    [~,density,xmesh,cdf] = kde(data(start_point(i):end_point(i)));
    [~,ind] = max(density);
    level(i) = xmesh(ind);
    cdf_val(i) = cdf(ind)*100;
end

% interp_points = start_point + round(ln_seg/2);
% 
% if length(level) > 1
%     baseline = interp1(interp_points,level,1:T,'spline');
% else
%     baseline = level*ones(1,T);
% end