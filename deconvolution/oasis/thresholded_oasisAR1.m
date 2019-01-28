function [c, s, b, g, smin, active_set] = thresholded_oasisAR1(y, g, sn, optimize_b,...
    optimize_g, decimate, maxIter, thresh_factor)
%% Infer the most likely discretized spike train underlying an AR(1) fluorescence trace
% Solves the sparse non-negative deconvolution problem
%  min 1/2|c-y|^2  subject to s_t = c_t-g c_{t-1} >=s_min or =0

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
%impulse response.
%   sn:  scalar, standard deviation of the noise distribution
%   optimize_b: bool, optimize baseline if True
%   optimize_g: integer, number of large, isolated events to consider for
%       optimizing g
%   decimate: int, decimation factor for estimating hyper-parameters faster
%       on decimated data
%   maxIter:  int, maximum number of iterations
%   active_set: npool x 4 matrix, warm stared active sets
%  thresh_factor: scalar, set the maximum thresh as thresh_factor*sn^2*T

%% outputs
%   c: T*1 vector, the inferred denoised fluorescence signal at each time-bin.
%   s: T*1 vector, discetized deconvolved neural activity (spikes)
%   b: scalar, fluorescence baseline
%   g: scalar, parameter of the AR(1) process
%   smin: scalar, minimum nonzero spike count 
%   active_set: npool x 4 matrix, active sets

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging


%% input arguments
y = reshape(y, [], 1);
T = length(y);

if ~exist('g', 'var') || isempty(g)
    g = estimate_time_constant(y, 1);
end
if ~exist('sn', 'var') || isempty(sn)
    sn = GetSn(y);
end
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
if ~exist('thresh_factor', 'var') || isempty(thresh_factor)
    thresh_factor = 1.0;
end

thresh = thresh_factor* sn * sn * T;
smin = 0;

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
if ~optimize_b   %% don't optimize the baseline b
    %% initialization
    b = 0;
    [solution, spks, active_set] = oasisAR1(y, g, [], smin);
            len_active_set = size(active_set, 1); 

    %% iteratively update parameters smin & g
    for miter=1:maxIter
        res = y - solution;
        RSS = res' * res;
        if or(RSS>thresh, sum(solution)<1e-9)  % constrained form has been found, stop
            break;
        else
            % update s_min
            [smin, solution, spks, active_set] = update_smin(y, g, smin,...
                solution, spks, active_set, sqrt(thresh)); 
      
            % update g
            if and(optimize_g, ~g_converged);
                g0 = g;
                [solution, active_set, g, spks] = update_g(y, active_set,smin);
                if abs(g-g0)/g0 < 1e-3 % g is converged
                    g_converged = true;
                end
            end
            
            % no more change of the active set 
            if size(active_set,1)==len_active_set
                break; 
            end
        end
    end
else
    %% initialization
    b = quantile(y, 0.15); 
    [solution, spks, active_set] = oasisAR1(y-b, g, [], smin);
    
    %% optimize the baseline b and dependends on the optimized g too
    g_converged = false;
    for miter=1:maxIter
        res = y - solution - b;
        RSS = res' * res;
        len_active_set = size(active_set,1);
        
        if or(abs(RSS-thresh) < tol, sum(solution)<1e-9)
            break;
        else
            % update smin
            [smin, solution, spks, active_set] = update_smin(y-b, g, smin,...
                solution, spks, active_set, sqrt(thresh)); 
            b = mean(y-solution); 
            
            % update b and g
            if and(optimize_g, ~g_converged);
                g0 = g;
                [solution, active_set, g, spks] = update_g(y-b, active_set, smin);       
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
    function [smin, solution, spks, active_set] = update_smin(y, g, smin, solution, ...
            spks, active_set, thresh) 
         %%estimate smin to match the thresholded RSS 
         len_active_set = size(active_set, 1); 
         s_max = max(active_set(:,1)./active_set(:,2)); 
         sv = linspace(smin, s_max, min(9, len_active_set));
         ind_start = 1;
         ind_end = length(sv);
         
         while (ind_end-ind_start)>1
             ind = floor((ind_start+ind_end)/2);
             tmp_smin = sv(ind);
             [tmp_solution, tmp_spks, tmp_active_set] = oasisAR1([], g, [], ...
                 tmp_smin, active_set);
             sqRSS = norm(y-tmp_solution,2);
             if sqRSS<thresh % increase smin
                 solution = tmp_solution; 
                 spks = tmp_spks; 
                 active_set = tmp_active_set; 
                 smin = tmp_smin; 
                 ind_start = ind;
             elseif sqRSS>thresh % decrease smin
                 ind_end = ind; 
             else
                 break; 
             end
         end
    end


end

 %update the AR coefficient: g
function [c, active_set, g, s] = update_g(y, active_set, smin)
%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities
%withone entry per time-bin.
%   active_set: npools*4 matrix, previous active sets.
% smin: scalr, minimize size of nonzero spikes 

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
yp = y; 
for m=1:len_active_set
    tmp_h = exp(log(g)*(0:maxl)');   % response kernel
    tmp_hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
    li = active_set(m, 4);
    ti = active_set(m, 3);
    idx = ti:(ti+li-1);
    active_set(m,1) = (yp(idx))'*tmp_h(1:li);
    active_set(m,2) = tmp_hh(li);
end
[c,s,active_set] = oasisAR1(y, g, 0, smin, active_set);

%% nested functions
    function rss = rss_g(g)
        h = exp(log(g)*(0:maxl)');   % response kernel
        hh = cumsum(h.*h);        % hh(k) = h(1:k)'*h(1:k)
        yp = y;     % include the penalty term
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