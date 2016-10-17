function [obj_ds, Y_ds] = downSample(obj, Y)
%% downsample data
% input:
%   Y: d*T matrix or d1*d2*T matrix
% output:
%   obj_ds: new object correspond to downsampled data
%   Y_ds:   d1s*d2s*Ts, downsampled data

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

Y_ds = dsData(Y, obj.options);    % call dsData to do downsampling
[d1, d2, ~] = size(Y_ds);
obj_ds = obj.copy();
obj_ds.Fs = obj.Fs/obj.options.tsub;
ssub = obj.options.ssub;
tsub = obj.options.tsub;
% update parameters
obj_ds.updateParams('d1', d1, 'd2',...
    d2, 'ssub', 1, 'tsub', 1, ...
    'gSig', ceil(obj.options.gSig/ssub), ...
    'gSiz', ceil(obj.options.gSiz/ssub), ...
    'bSiz', ceil(obj.options.bSiz/ssub));

obj_ds.options.Fs = obj.Fs/tsub;
if tsub~=1
    obj_ds.kernel = dsKernel(obj.kernel, tsub); 
end 

end