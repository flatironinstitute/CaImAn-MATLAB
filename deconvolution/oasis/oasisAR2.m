function [c, s, active_set] = oasisAR2(y, g, lam, smin, T_over_ISI, jitter, active_set)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g1*c_{t-1}-g2*c_{t-2} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g1:  scalar, first parameter of the AR(2) process that models the fluorescence ...
%impulse response.
%   g2:  scalar, second parameter of the AR(2) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, sparsity penalty parameter lambda.
%   smin: scalar, optional, default 0
%miniumal non-zero activity within each bin (minimal 'spike size').
%   T_over_ISI: scalar, ratio of recording duration T and maximumal
%       inter-spike-interval. default: 1
%   jitter: bool, perform correction step by jittering spike times to
%       minimize RSS. it helps to avoid delayed spike detection. default:
%       false;

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   active_set: npool * 4 matrix, active set

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
y = reshape(y, [], 1);
if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y, 2);
end
g1 = g(1); 
g2 = g(2); 

if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('smin', 'var') || isempty(smin);   smin = 0; end
if ~exist('T_over_ISI', 'var') || isempty(T_over_ISI)
    T_over_ISI = 1;
end
if ~exist('jitter', 'var') || isempty(jitter)
    jitter = false;
end

%% initialization
T = length(y);
yp = y - lam * (1-g1-g2);
yp(end-1) = y(end-1) - lam*(1-g1);
yp(end) = y(end) - lam;

if ~exist('active_set', 'var') || isempty(active_set)
    % active set
    len_active_set = length(y);
    active_set = [yp, yp, (1:T)', ones(T,1), (1:T)'-1, (1:T)'+1];
    active_set(1,5) = nan;
    active_set(end,6) = nan;
else
    len_active_set = size(active_set,1);
    active_set(:,5) = [nan;(1:len_active_set-1)'];
    active_set(:,6) = [(2:len_active_set)';nan];
end
idx = true(len_active_set,1);

% precompute
len_g = T / T_over_ISI;
temp = roots([1, -g1, -g2]);
d = max(temp); r = min(temp);
g11 = (exp(log(d)*(1:len_g)') - exp(log(r)*(1:len_g)')) / (d-r);
g12 = [0; g2*g11(1:(end-1))];
g11g11 = cumsum(g11.*g11);
g11g12 = cumsum(g11.*g12);

%% run OASIS
ii = 2;
ii_next = active_set(ii,6);
ii_prev = active_set(ii,5);
while ~isnan(ii_next)
    % find the active set
    while (~isnan(ii_next)) && (g11(active_set(ii,4))*active_set(ii,1) + ...
            g12(active_set(ii,4))*active_set(ii_prev,2) + smin <= active_set(ii_next,1))
        active_set(ii_next,5) = ii;         
        ii = ii_next;       % move to the next pool 
        ii_next = active_set(ii,6);     
        ii_prev = active_set(ii,5);
    end
    if isnan(ii_next); break; end
    
    %% merge pools
    active_set(ii, 4) = active_set(ii, 4) + active_set(ii_next, 4);
    ti = active_set(ii,3);
    li = active_set(ii,4);
    active_set(ii,1) = (g11(1:li)'*yp(ti+(1:li)-1) - g11g12(li)*active_set(ii_prev,2))/g11g11(li);
    active_set(ii,2) = g11(li)*active_set(ii,1) + g12(li)*active_set(ii_prev,2);
    idx(ii_next) = false; 
    active_set(ii,6) = active_set(ii_next, 6); 
    ii_next = active_set(ii, 6); 
    if ~isnan(ii_next)
        active_set(ii_next,5) = ii; 
    end
    
    ii_prev_1 = active_set(ii,5);
    ii_prev_2 = active_set(ii_prev_1,5);
    %% backtrack until violations fixed
    while (~isnan(ii_prev_2)) &&  (g11(active_set(ii_prev_1,4))*active_set(ii_prev_1,1) + g12(active_set(ii_prev_1,4))...
            *active_set(ii_prev_2,2) + smin > active_set(ii,1))
        ii_next = ii; 
        ii = ii_prev_1; 
        ii_prev = ii_prev_2; 
        
        active_set(ii, 4) = active_set(ii, 4) + active_set(ii_next, 4);
        ti = active_set(ii,3);
        li = active_set(ii,4);
        active_set(ii,1) = (g11(1:li)'*yp(ti+(1:li)-1) - g11g12(li)*active_set(ii_prev,2))/g11g11(li);
        active_set(ii,2) = g11(li)*active_set(ii,1) + g12(li)*active_set(ii_prev,2);
        idx(ii_next) = false;
        active_set(ii,6) = active_set(ii_next, 6);
        ii_next = active_set(ii, 6);
        if ~isnan(ii_next)
            active_set(ii_next,5) = ii;
        end
        ii_prev_1 = active_set(ii,5);
        ii_prev_2 = active_set(ii_prev_1,5);
    end
    
end

%% jitter
% a_s = active_set;
if jitter
    disp('to be done\n');
end

active_set(~idx, :) = [];
len_active_set = size(active_set,1);

%% construct solution for all t
c = zeros(T,1);
for ii=1:len_active_set
    vi = active_set(ii,1);
    ti = active_set(ii,3);
    li = active_set(ii,4);
    c(ti) = vi;
    
    for m=1:(li-1)
        c(ti+m) = g1*c(ti+m-1) + g2*c(ti+m-2);
    end
end
c(c<0) = 0;

s = [0;0;0; c(4:end)-g1*c(3:(end-1))-g2*c(2:(end-2))];
s(s<smin) = 0;
