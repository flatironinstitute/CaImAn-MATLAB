function [c, s, b, g, lam, active_set] = constrained_oasisAR2(y, g, sn, optimize_b,...
    optimize_g, decimate, maxIter)
%% Infer the most likely discretized spike train underlying an AR(2) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min lam |s|_1
%  subject to |y-c|_2^2 <= sn^2*T

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  2 x 1 vector, Parameter of the AR(2) process that models the fluorescence ...
%impulse response.
%   sn:  scalar, standard deviation of the noise distribution
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data
%   maxIter:  int, maximum number of iterations
%   active_set: npool x 4 matrix, warm stared active sets

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   lam: scalar, sparsity penalty parameter
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging


%% input arguments
y = reshape(y, [], 1);
T = length(y);

if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y, 2);
end
if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(y);
end
if ~exist('lam', 'var') || isempty(lam);   lam = 0; end
if ~exist('optimize_b', 'var') || isempty(optimize_b)
    optimize_b = false;
end
if ~exist('optimize_g', 'var') || isempty(optimize_g)
    optimize_g = 0;
end
if ~exist('decimate', 'var') || isempty(decimate)
    decimate = 1;
else
    decimate = max(1, round(decimate));
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 10;
end

thresh = sn * sn * T;
lam = 0;

% change parameters due to downsampling
if decimate>1
    decimate = 1;  %#ok<NASGU>
    disp('to be done');
    %     fluo = y;
    %     y = resample(y, 1, decimate);
    %     g = g^decimate;
    %     thresh = thresh / decimate / decimate;
    %     T = length(y);
end
g_converged = false;

%% optimize parameters
tol = 1e-4;
flag_lam = true;
if ~optimize_b   %% don't optimize the baseline b
    %% initialization
    b = 0;
    [solution, spks, active_set] = oasisAR2(y, g, lam);
    
    %% iteratively update parameters lambda & g
    for miter=1:maxIter
        res = y - solution;
        RSS = res' * res;
        len_active_set = size(active_set, 1); 
        if or(abs(RSS-thresh) < tol, ~flag_lam)  % constrained form has been found, stop
            break;
        else
            % update lam
            update_phi; 
            lam = lam + dphi; 
      
            % update g
            if and(optimize_g, ~g_converged);
                g0 = g;
                [solution, active_set, g, spks] = update_g(y, active_set,lam);
                if abs(g-g0)/g0 < 1e-3  % g is converged
                    g_converged = true;
                end
            end
        end
    end
else
    %% initialization
    b = quantile(y, 0.15); 
    [solution, spks, active_set] = oasisAR2(y-b, g, lam);
    update_lam_b; 
    
    %% optimize the baseline b and dependends on the optimized g too
    g_converged = false;
    for miter=1:maxIter
        res = y - solution - b;
        RSS = res' * res;
        len_active_set = size(active_set,1);
        
        if or(abs(RSS-thresh) < tol, sum(solution)<1e-9)
            break;
        else
            %% update b & lamba
            update_phi();
            update_lam_b();
                        % update b and g    
            % update b and g
            if and(optimize_g, ~g_converged);
                g0 = g;
                [solution, active_set, g, spks] = update_g(y-b, active_set,lam);       
                if abs(g-g0)/g0 < 1e-4;
                    g_converged = true;
                end
            end

        end
    end
    
end
c = solution;
s = spks;

%% nested functions 
    function update_phi()  % estimate dphi to match the thresholded RSS
        zeta = zeros(size(solution));
        maxl = max(active_set(:, 4));
        h = g.^(0:maxl);
        for ii=1:len_active_set
            ti = active_set(ii, 3);
            li = active_set(ii, 4);
            idx = 0:(li-1);
            if ii<len_active_set
                zeta(ti+idx) = (1-g^li)/ active_set(ii,2) * h(1:li);
            else
                zeta(ti+idx) = 1/active_set(ii,2) * h(1:li);
            end
        end
        
        if optimize_b
            zeta = zeta - mean(zeta);
            tmp_res = res - mean(res);
            aa = zeta' * zeta;
            bb = tmp_res' * zeta;
            cc = tmp_res'*tmp_res - thresh;
            dphi = (-bb + sqrt(bb^2-aa*cc)) / aa;
        else
            aa = zeta'*zeta;
            bb = res'*zeta;
            cc = RSS-thresh;
            dphi = (-bb + sqrt(bb^2-aa*cc)) / aa;
        end
        if imag(dphi)>1e-9
            flag_phi = false; 
            return; 
        else
            flag_phi = true; 
        end
        active_set(:,1) = active_set(:,1) - dphi*(1-g.^active_set(:,4));
        [solution, spks, active_set] = oasisAR2([], g, lam, [], active_set);
    end

    function update_lam_b() % estimate lambda  & b
        db = mean(y-solution) - b;
        b = b + db;
        dlam = -db/(1-g);
        
        lam = lam + dlam;
        % correct the last pool
        active_set(end,1) = active_set(end,1) - lam*g^(active_set(end,4));
        ti = active_set(end,3); li = active_set(end,4); idx = 0:(li-1);
        solution(ti+idx) = max(0, active_set(end,1)/active_set(end,2)) * (g.^idx);
    end

end
 %update the AR coefficient: g
function [c, active_set, g, s] = update_g(y, active_set, lam)
%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   active_set: npools*4 matrix, previous active sets.
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
[c,s,active_set] = oasisAR2(y, g, lam, [], active_set);

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