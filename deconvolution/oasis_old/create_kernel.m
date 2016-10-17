function kernel = create_kernel(kernel_type, pars, nMax, lb, ub, bound_pars)
%% create convolution kernel
%% inputs:
%   kernel_type: string, convolution kernel type. now support {'exp',
%       'exp2', 'vector'}
%   pars: parameters for the selected kernel type
%   nMax: length of the kernel
%   lb:     lower bound for each parameter
%   ub:     upper bound for each parameter
%   bound_pars: logical variable, bound the parameters or not {1, 0}
%% outputs
%   kernel: struct variable

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

%% kernel size 
if ~exist('nMax', 'var') || isempty(nMax)
    nMax = 50;
end
kernel.nMax = nMax; 

%% initialize kernel 
if ~exist('kernel_type', 'var') || isempty(kernel_type)
    kernel_type = 'exp2'; 
end
if strcmpi(kernel_type, 'exp')
    % single exponential function: ~ exp(-t/tau)
    kernel.type = 'exp';
    % parameter
    if ~exist('pars', 'var') || isempty(pars)
        kernel.pars = [5, .1];
    else
        kernel = kernel.pars; 
    end    
    % function handle 
    kernel.fhandle = @(pars, t) exp(-t/pars) * (1-exp(-1/pars)) ...
        / (1-exp(-nMax/pars)); 
elseif strcmpi(kernel_type, 'vector')
    % single vector 
    kernel.type = 'vector'; 
    % parameter 
    if ~exist('pars', 'var') || isempty(pars)
        kernel.pars = exp(-(1:nMax)/10); 
    else
        kernel.pars = pars; 
    end
    % function handle 
    kernel.fhandle = @(pars, t) pars/sum(pars); 
else
    % differencing of two exponential function:
    % ~  exp(-t/tau_d)-exp(-t/tau_r)
    kernel.type = 'exp2';
    
    % parameters
    if ~exist('pars', 'var') || isempty(pars)
        kernel.pars = [10, 1];
    else
        kernel.pars = pars; 
    end
    % function handle 
    kernel.fhandle = @(pars, t) (exp(-t/pars(1)) - exp(-t/pars(2)))  ...
        /( (1-exp(-nMax/pars(1)))/(1-exp(-1/pars(1))) ...
        - (1-exp(-nMax/pars(2)))/(1-exp(-1/pars(2))));
end

% lower and upper bounds for parameters 
if ~exist('lb', 'var') || isempty(lb)
    kernel.lb = 0.5*kernel.pars;
else
    kernel.lb = lb; 
end
if ~exist('ub', 'var') || isempty(ub)
    kernel.ub = 2*kernel.pars;
else
    kernel.ub = ub; 
end

% bound the parameters of not 
if ~exist('bound_pars', 'var')||isempty(bound_pars)
    kernel.bound_pars = false; 
else
    kernel.bound_pars = bound_pars; 
end