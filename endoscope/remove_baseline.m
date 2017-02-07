function [y, b] = remove_baseline(y, sn)
% estiamte baseline of the calcium traces and subtract it 
if ~exist('sn', 'var') || isempty(sn)
    sn = get_noise_fft(reshape(y, 1, [])); 
end
sz = size(y); 
y = reshape(y, 1, []); 
y_diff = [-1, diff(y)]; 
b = median(y(and(y_diff>=0, y_diff<sn))); 
y = reshape(y-b, sz); 