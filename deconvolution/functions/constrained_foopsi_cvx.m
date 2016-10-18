function [c, s] = constrained_foopsi_cvx(y, g, sn, solver)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min |s|_1 subject to s=G*c >=0 and |y-c|_2^2=\sigma*T^2 

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   sn:  scalar, noise standard deviation 
%   solver: string, optimization solver

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% initialization
y = reshape(y, [], 1);
T = length(y);
if ~exist('g', 'var') || isempty(g);   g = 0.95; end
if ~exist('sn', 'var') || isempty(sn);   sn = GetSn(y); end
if ~exist('solver', 'var') || isempty(smin);   solver = 'SDPT3'; end

% construct the deconvolution matrix (s=G*c)
len_g = length(g);
g = reshape(g, 1, []);
indc = bsxfun(@minus, (1:T)', 0:len_g);
indr = repmat((1:T)', [1, len_g+1]);
v = ones(T,1)*[1,-g];
ind = (indc>0);
G = sparse(indr(ind), indc(ind), v(ind), T, T);

%% run optimization
cvx_solver(solver); 
cvx_begin quiet
    variable c(T)
    minimize(norm(c,1))
    subject to
        G*c >=0;
        norm(y-c,2) <= sn*sqrt(T); 
cvx_end

s = reshape(G*c, 1, []);
s(1) = 0;
c = reshape(c,1, []);



