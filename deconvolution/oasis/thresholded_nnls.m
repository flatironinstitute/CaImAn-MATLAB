function [c, s] = thresholded_nnls(y, g, sn, smin, shift, win, tol, maxIter, mask, threshold_factor)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|Ks-y|^2 + lam * |s|_1 subject to s_t = c_t-g c_{t-1} >= 0

%% inputs:
%   y:  T*1 vector, vth dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
%   g:  vector, shape (p,)
%       if p in (1,2): AR coefficients for AR(p) process 
%       else: kernel that models the fluorescence implulse response 
%   sn:  scalar, standard deviation of the noise 
%   smin: scalar, minimum spike size 
%   shift: integer scalar, number of frames by which to shift window from on run of
%       NNLS to the next, default-100
%   win: integer acalar, window size 
%   tol: scalar, tolerance parameters 
%   maxIter: scalar, maximum number of iterations before termination 
%   mask: T * 1 boolean vector, restrict potential spike times 
%   smin: scalar, minimum spike size 
%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes) 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% input arguments  
T = length(y); 
y = reshape(y, [], 1); 

if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(y);
end
if ~exist('shift', 'var') || isempty(shift)
    shift = 100; 
end
if ~exist('win', 'var') || isempty(win)
    win = 200; 
end

if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9; 
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = []; 
end
if ~exist('mask', 'var') || isempty(mask)
    mask = true(T,1); 
end
if ~exist('smin', 'var') || isempty(smin)
    smin = 0; 
end

thresh = thresh_factor* sn * sn * T;

%% get the response kernel
w = win; 
K = zeros(w); 
[u, t] = meshgrid(1:w, 1:w);  
ind = 1+t-u;
if length(g)==1
    h = exp(log(g)*(0:(w-1)));
elseif length(g)==2
    temp = roots([1, -g(1), -g(2)]);
    d = max(temp);
    r = min(temp);
    h = (exp(log(d)*(1:w)) - exp(log(r)*(1:w))) / (d-r); % convolution kernel
else
    h = g;
end
K(ind>0) = h(ind(ind>0));   % convolution matrix
KK = K'*K;

%% initialization 
a = sum(inv(K));
yp = y - lam * a(1);
yp((end-w+1):end) = yp((end-w+1):end) - lam * a';
s = zeros(T, 1); 
c = zeros(T, 1); 

%% run online deconvolution 
t = 1; 
yp0 = yp; 
while t <= T-w+1
    ind = t:(t+w-1); 
    s(ind) = nnls(KK, K'*yp(ind), s(ind), tol, maxIter, mask(ind)); 
    yp(ind) = yp(ind) - K(:, 1:shift)*s(t:(t+shift-1)); 
    c(ind) = c(ind) + K(:, 1:shift)*s(t:(t+shift-1)); 
    t = t + shift; 
end 
s(t:T) = nnls(KK((t+w-T):w, (t+w-T):w), K(1:(T-t+1), 1:(T-t+1))'*yp(t:T), ...
    s(t:T), tol, maxIter, mask(t:T)); 
c(t:T) = c(t:T) + K((t+w-T):w, (t+w-T):w) * s(t:T); 

%% running thresholded version for fast initialization 
if smin>0
    yp = yp0; 
    c = zeros(size(c)); 
    mask = (s>1e-4); 
    t = 1;
    while true
        [tmp_min, ind] = min(s(mask)); 
        ind = t:(t+w-1);
        s(ind) = nnls(KK, K'*yp(ind), s(ind), tol, maxIter, mask(ind), smin);
        yp(ind) = yp(ind) - K(:, 1:shift)*s(t:(t+shift-1));
        c(ind) = c(ind) + K(:, 1:shift)*s(t:(t+shift-1));
        t = t + shift;
    end
    s(t:T) = nnls(KK((t+w-T):w, (t+w-T):w), K(1:(T-t+1), 1:(T-t+1))'*yp(t:T), ...
        s(t:T), tol, maxIter, mask(t:T), smin);
    c(t:T) = c(t:T) + K((t+w-T):w, (t+w-T):w) * s(t:T); 
end

%% 

function s = nnls(KK, Ky, s, tol, maxIter, mask, smin)
%% fast algorithm for solving nonnegativity constrained least squared
% problem minize norm(y-K*s, 2), s.t. s>=0. 

%% inputs: 
%   KK: p x p matrix, K'*K
%   Ky: n x 1 vector, K'*y
%   s: p x 1 vector, warm started s 
%   tol: scalar, smallest nonzero values 
%   maxIter: scalar, maximum nonzero values 
%   mask: p x 1 vector, mask to restrict potential spike times considered
%   smin: scala, minimize size of the spike 

%% outputs: 
%   s: p x 1 vector, solution 

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging
% Bro R & Jong S, Journal of Chemometrics 1997, A FAST NON-NEGATIVITY-CONSTRAINED LEAST SQUARES ALGORITHM


%% input arguments 
if ~exist('mask', 'var')||isempty(mask)
    mask = true(size(KK,1), 1); 
elseif any(mask)
    KK = KK(mask, mask); 
    Ky = Ky(mask); 
else
    s = double(mask); 
    return; 
end

p = size(KK,2);      % dimension of s 
if ~exist('smin', 'var') || isempty(smin)
    smin = 0; 
end
vth = ones(p,1)*smin;  
if ~exist('s', 'var') || isempty(s)
    s = zeros(p, 1); 
    l = Ky; 
    Pset = false(p,1); 
else
    s = s(mask); 
    s = s - smin; 
    Pset = (s>0);    
    s(~Pset) = 0; 
    l = Ky - KK*s - KK * Pset*smin; 
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9; 
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = p; 
end


% outer loop: loop for passive set
[~, ind] = max(l); 
for miter =1:maxIter   
    % remove element from the active set 
    Pset(ind) = true;
    
    % solve unconstrained least squares over the passive set 
    try
        mu = KK(Pset, Pset) \ (Ky(Pset)-KK(Pset, Pset)*vth(Pset));
    catch % sigular issue
        mu = (KK(Pset, Pset) + tol*eye(sum(Pset))) \ (Ky(Pset)...
            -KK(Pset, Pset)*vth(Pset));
    end
    
    % inner loop: correct nonnegativity violations 
    while any(mu<=0)   
        temp = s(Pset) ./ (s(Pset)-mu);
        a = min(temp(mu<0)); 
        s(Pset) = s(Pset) + a*(mu-s(Pset));
        Pset(s<=tol) = false;
        
        % solve unconstrained least squares over the passive set
        try
            mu = KK(Pset, Pset) \ (Ky(Pset)-KK(Pset, Pset)*vth(Pset));
        catch % sigular issue
            mu = (KK(Pset, Pset) + tol*eye(sum(Pset))) \ (Ky(Pset)...
                -KK(Pset, Pset)*vth(Pset));
        end
    end
    
    s(Pset) = mu; 
    l = Ky - KK*s - KK*Pset*smin; 
        % at least one iteration 
    [lmax, ind] = max(l); 
    if lmax < tol         % no more passive set
        break;
    end
end
s(Pset) = mu + smin; 
temp = double(mask);
temp(mask) = s; 
s = temp; 
