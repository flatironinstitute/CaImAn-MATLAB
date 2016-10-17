function [c, active_set, g, s] = update_g(y, active_set, g, lam)
%% update the tuning parameter lambda to reduce |s|_1 while |y-c|_2^2 <= sn^2*T

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   active_set: npools*4 matrix, previous active sets
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   lam:  scalar, curret value of sparsity penalty parameter lambda.

%% outputs
%   c: T*1 vector
%   s: T*1 vector, spike train
%   active_set: npool x 4 matrix, active sets
%   g: scalar

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

%% initialization
len_active_set = size(active_set, 1);  %number of active sets
y = reshape(y,[],1);    % fluorescence data
maxl = max(active_set(:, 4));   % maximum ISI
c = zeros(size(y));     % the optimal denoised trace

%% find the optimal g and get the warm started active_set
g = fminbnd(@rss_g, 0, 1);
yp = y - lam*(1-g); 
for m=1:len_active_set
    tmp_h = exp(log(g)*(0:maxl)');   % response kernel
    tmp_hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
    li = active_set(m, 4);
    ti = active_set(m, 3);
    idx = ti:(ti+li-1);
    active_set(m,1) = (yp(idx))'*tmp_h(1:li);
    active_set(m,2) = tmp_hh(li);
end
[c,s,active_set] = oasisAR1(y, g, lam, [], active_set);

%% nested functions
    function rss = rss_g(g)
        h = exp(log(g)*(0:maxl)');   % response kernel
        hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
        yp = y - lam*(1-g);     % include the penalty term
        for ii=1:len_active_set
            li = active_set(ii, 4);
            ti = active_set(ii, 3);
            idx = ti:(ti+li-1);
            tmp_v = max(yp(idx)' * h(1:li) / hh(li), 0);
            c(idx) = tmp_v*h(1:li);
        end
        res = y-c;
        rss = res'*res;     % residual sum of squares
    end
end
























