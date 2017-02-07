function kernel_new = dsKernel(kernel, tsub)
%% downsample/upsample the convolution kernel 
%% inputs:
%   kernel: struct variable with fields {'kernel_type', 'pars', 'nMax', 'lb', 'ub', 'bound_pars'}
%       kernel_type: string, convolution kernel type. now support {'exp',
%       'exp2', 'vector'}
%       pars: parameters for the selected kernel type
%       nMax: length of the kernel
%       lb:     lower bound for each parameter
%       ub:     upper bound for each parameter
%       bound_pars: logical variable, bound the parameters or not {1, 0}
%% outputs
%   kernel: struct variable

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

kernel_new = kernel; 
if nargin<2 || isempty(tsub)
    return;
end

%% kernel size 
kernel_new.nMax = ceil(kernel.nMax/tsub); 

%% kernel type 
kernel_type = kernel.type; 
if strcmpi(kernel_type, 'exp')
    % single exponential function: ~ exp(-t/tau)
    kernel_new.pars = kernel.pars/tsub; 
elseif strcmpi(kernel_type, 'vector')
    % single vector 
len_kernel = length(kernel.pars); 
    kernel_new.pars = resample(kernel.pars, ceil(len_kernel/tsub), len_kernel); 
else
    % differencing of two exponential function:
    % ~  exp(-t/tau_d)-exp(-t/tau_r)
    kernel_new.pars = kernel.pars/tsub; 
end

%% lower and upper bounds for parameters 
kernel_new.lb = kernel.lb / tsub; 
kernel_new.ub = kernel.ub / tsub; 