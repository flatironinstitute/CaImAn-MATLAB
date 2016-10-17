function [c, s, active_set] = oasisAR1(y, g, lam, smin, active_set)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
% OR %%
% len_active_set*4 matrix, active set

%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, sparsity penalty parameter lambda.
%   smin: scalar, optional, default 0
%miniumal non-zero activity within each bin (minimal 'spike size').
%   active_set: npool x 4 matrix, warm stared active sets

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
y = reshape(y, [], 1);
T = length(y);
if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y);
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('smin', 'var') || isempty(smin);   smin = 0; end
if ~exist('active_set', 'var') || isempty(active_set)
    len_active_set = T;
    active_set = [y-lam*(1-g),ones(T,1),(1:T)',ones(T,1),(1:T)'-1, (1:T)'+1]; % each row is one pool: (vi, wi, t, l)
    active_set(end, :) = [y(end)-lam,1,T,1,T-1,nan] ;
    active_set(1,5) = nan;
else
    len_active_set = size(active_set,1);
    active_set(:,5) = [nan; (1:len_active_set-1)']; 
    active_set(:,6) = [(2:len_active_set)';nan]; 
end
idx = true(len_active_set,1);

%% run OASIS
ii = 1;
ii_next = active_set(ii,6);
while ~isnan(ii_next)
    % find the active set
    while (~isnan(ii_next)) && (active_set(ii_next,1)/active_set(ii_next,2)...
            >=active_set(ii,1)/active_set(ii,2)*g^(active_set(ii,4))+smin)
        active_set(ii_next,5) = ii;
        ii = ii_next; 
        ii_next = active_set(ii,6);
    end
    
    if isnan(ii_next); break; end
    
    %% merge pools
    active_set(ii,1) = active_set(ii,1) + active_set(ii_next,1)* (g^(active_set(ii,4)));
    active_set(ii,2) = active_set(ii,2) + active_set(ii_next,2)*(g^(2*active_set(ii,4)));
    active_set(ii,4) = active_set(ii,4) + active_set(ii_next,4);
    active_set(ii,6) = active_set(ii_next,6);
    idx(ii_next) = false;
    ii_next = active_set(ii,6);
    ii_prev = active_set(ii, 5);

    %% backtrack until violations fixed
    while (~isnan(ii_prev)) && (active_set(ii,1)/active_set(ii,2)<...
            active_set(ii_prev,1)/active_set(ii_prev,2)*g^(active_set(ii_prev,4))+smin)
        ii_next = ii;
        ii = ii_prev;
        active_set(ii,1) = active_set(ii,1) + active_set(ii_next,1)* (g^(active_set(ii,4)));
        active_set(ii,2) = active_set(ii,2) + active_set(ii_next,2)*(g^(2*active_set(ii,4)));
        active_set(ii,4) = active_set(ii,4) + active_set(ii_next,4);
        active_set(ii,6) = active_set(ii_next,6);
        idx(ii_next) = false;
        
        ii_prev = active_set(ii, 5);
        ii_next = active_set(ii,6);
    end
end
active_set(~idx, :) = [];
len_active_set = size(active_set,1);

%% construct solution for all t
c = zeros(size(y));
s = c;
for ii=1:len_active_set
    t0 = active_set(ii,3);
    tau = active_set(ii, 4);
    c(t0:(t0+tau-1)) = max(0,active_set(ii,1)/active_set(ii,2)) * (g.^(0:(tau-1)));
end

s(active_set(2:end,3)) = c(active_set(2:end,3)) - g*c(active_set(2:end,3)-1);
