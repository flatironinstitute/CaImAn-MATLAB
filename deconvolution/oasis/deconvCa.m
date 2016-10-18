function [c, s, kernel, iter] = deconvCa(y, kernel, smin, fit_gt, debug_on, sn, maxIter, theta, lambda)
%% deconvolve calcium traces to infer spike counts
%% inputs:
%   y:  1*T vector, observed calcium traces
%   kernel: struct variable with two fields {'fhandle', 'pars', nMax}. kernel
%   deterines the convolution kernel
%      fhandle: function handle, the function form of the response function
%       pars: p*1 vector, parameters for fhandle
%       nMax: scalar, maximum number of frames required by calcium indicators
%           to return resting states
%   smin: scalar, minimum number of nonzero spike count within each bin
%   fit_gt: true or false. iteratively fit gt or not
%   debug_on: play and save video of the whole procedure.
%   sn:     noise power
%   maxIter: scalar, maximum iterations for running OASIS, default (3)
%   theta: 1*T vecotr,  weight vector at each time point
%   lambda: scalar, tuning parameter

%% outputs:
%   c: 1*T vector, inferred calcium trace
%   s: 1*T vector, inferred spike train
%   kernel: convolution kernel
%   iters: number of iterations for updating irf

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
% This work is based on one NIPS paper by Johannes Friedrich & Liam
% Paninski

%% input arugments
y = reshape(y, 1, []);  % convert data into row vector
if ~exist('sn', 'var') || isempty(sn)
    sn = get_noise_fft(y);
end
T = length(y);          % number of frames

% get convolution kernel
if ~exist('kernel', 'var') || isempty(kernel)
    kernel = create_kernel('exp2');
end
fhandle = kernel.fhandle;
nMax = kernel.nMax;

y = [y, zeros(1, nMax)];  % add few more elements for computation convenience

% threshold of the spike
if ~exist('smin', 'var') || isempty(smin)
    smin = 3*sn;               % min spike count
else
    smin = smin * sn;
end

% fit response function or not
if ~exist('fit_gt', 'var') || isempty(fit_gt)
    fit_gt = false;
end
thresh = 1e-3;

% debug mode 
if ~exist('debug_on', 'var') || isempty(debug_on)
    debug_on = false;
end

% maximum iterations for running OASIS 
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 5; 
end

% tuning parameter for enforcing sparsity
if ~exist('lambda', 'var') || isempty(lambda)
    lambda = 0;   % tuning parameter
end

% weight for each frame
if ~exist('theta', 'var')|| isempty(theta) || (length(theta)==1)     % weight vector
    theta = ones(1, T);
elseif length(theta)==T
    theta = reshape(theta, 1, T);
else
    disp('elements in theta should be equal to elements in y');
    return;
end

% tuning parameter for enforcing sparsity
if ~exist('lambda', 'var') || isempty(lambda)
    lambda = 0;   % tuning parameter
end

%% running OASIS
if debug_on
    figure('position', [1, 1 1500, 250]);
    plot(y); hold on;
    a = plot(y, 'r');
    avi_file = VideoWriter('example.avi');
    avi_file.open();
    xlim([1, T]);
end

f0 = inf;       % objective function
y0 = y;         % backup the raw trace
iter = 1;       % iteratiosn
for iter = 1:maxIter
    if debug_on
        title(sprintf('Iter %d, norm(residual) = %.4f', iter, f0));
    end
    % normalize gt to make its maximum value to be 1
    gt = fhandle(kernel.pars, 1:nMax); % 1*nMax, response function of calcium transients
    ind0 = 1; %correct the response function of the first event
    gt = [reshape(gt, 1, nMax), zeros(1, T)]; % add few more elements to gt for computation convenience
    
    gt_max = max(gt);
    gt = gt/gt_max;
    gt1 = [gt(ind0:end), zeros(1, ind0-1)]; % kernel for the first frame
    
    %% initialize values for the results
    s = zeros(size(y));     % spike count within each bin
    v = zeros(size(y));     % sum(theta_t*y_t) - lambda
    w = zeros(size(y));     % sum(theta_t^2 * g(t-ti+1)^2))
    pool_ti = zeros(size(y)); % first frame of each pool
    
    % initialize the first pool
    v(1) = theta(1)*y(1)*gt1(1) - lambda;
    w(1) = (theta(1)*gt1(1))^2;
    s(1) = max(0, v(1) / w(1));
    %     if s(1)<smin
    %         s(1) = 0;
    %     end
    pool_ti(1) = 1;
    tpre = 1;   % time of the last event
    
    %% start to initialize the next pool
    ii = 2;     % start to initialize the next pool
    t = 2;      % frame to begin
    frame_valid = true;
    while t<=T
        % plot results
        if debug_on
            c = zeros(1,T);
            c(1:T) = s(1) * gt1(1:T);
            for m=2:(ii-1)
                t0 = pool_ti(m);
                c(t0:T) = c(t0:T)+s(t0) * gt(1:(T-t0+1));
            end
            delete(a);
            a = plot(c(1:t), 'r', 'linewidth', 2);
            drawnow;
            temp = getframe();
            temp.cdata = imresize(temp.cdata, [200, 1500]);
            avi_file.writeVideo(temp);
        end
        %
        % find the spike count that minimizes the objective function
        if (t-tpre)>=nMax  % no influences from the previous events
            vi = theta(t)*y(t)*gt(1) - lambda;
        elseif tpre>1        % estimate vi and subtract the effect of the previous event
            vi = theta(t) * (y(t)-s(tpre)*gt(t-tpre+1)) * gt(1) - lambda;
        else   % special treatment with the first event
            vi = theta(t) * (y(t)-s(tpre)*gt1(t-tpre+1)) * gt(1) - lambda;
        end
        wi = (theta(t) * gt(1))^2;
        si = vi/wi;
        
        % check the violatoin of si
        if si>smin && frame_valid % no violation, create a new pool to save the result and move to the next new pool
            % peel off the previous event
            if (t-tpre)<nMax
                if tpre>1
                    y(tpre+(1:nMax)-1) = y(tpre+(1:nMax)-1) - s(tpre)*gt(1:nMax);
                else
                    y(tpre+(1:nMax)-1) = y(tpre+(1:nMax)-1) - s(tpre)*gt1(1:nMax);
                end
            end
            pool_ti(ii) = t;   % create a new pool
            tpre = t;
            v(t) = vi;
            w(t) = wi;
            s(t) = si;
            ii = ii+1;  % move to the next pool
            t = t+1;
        else  % with violation
            frame_valid = true; % allows the next frame to be valid
            if (t-tpre)>=nMax
                % ignore this frame directly and move to the next frame
                t = t + 1;
                continue;
            elseif tpre>1         % merge it to the current pool
                v(tpre) = v(tpre) + theta(t)*y(t)*gt(t-tpre+1);
                w(tpre) = w(tpre) + (theta(t)*gt(t-tpre+1))^2;
            else
                v(tpre) = v(tpre) + theta(t)*y(t)*gt1(t-tpre+1);
                w(tpre) = w(tpre) + (theta(t)*gt1(t-tpre+1));
            end
            s(tpre) = v(tpre)/w(tpre);  % update the current event
            
            % check the violation of the current pool
            if s(tpre)>smin      % the previous event is still avaiable
                t = t+1;
                continue;
            elseif ii==2  %the previous pool is not available anymore, but it's the first pool
                s(tpre) = max(0, s(tpre));
                t = t+1;
            else   % not available, then delete the current pool and force its first frame to be
                % event-free.
                t = tpre;  % go back to the first frame of the current pull
                
                frame_valid = false; % force the frame t to be invalid
                ii = ii-1;
                tpre = pool_ti(ii-1);
                
                % add back the signal of the previous pool
                if (t-tpre)>nMax
                    t = t+1;
                    frame_valid = true;
                    continue;
                elseif tpre>1
                    y(tpre+(1:nMax)-1) = y(tpre+(1:nMax)-1) + s(tpre)*gt(1:nMax);
                else
                    y(tpre+(1:nMax)-1) = y(tpre+(1:nMax)-1) + s(tpre)*gt1(1:nMax);
                end
            end
        end
    end
    %% collect the results
    pool_ti(ii:end) = [];
    temp = s(pool_ti);
    s = zeros(1, T);
    s(pool_ti) = temp;
    temp = s;
    temp(1) = 0;
    c = conv(temp, gt(1:nMax));
    c = c(1:T) + s(1)*gt1(1:T);
    f1 = norm(y(1:T)-c, 2);
    s = s*sum(gt);
    if debug_on
        delete(a);
        a = plot(c, 'r');
        drawnow();
        temp = getframe();
        temp.cdata = imresize(temp.cdata, [200, 1500]);
        avi_file.writeVideo(temp);
    end
    
    %% break the while loop
    if ~fit_gt  % don't iterate gt
        break;
    elseif (f0-f1)/f0 <= thresh % improvement is small, stop
        break;
    else  % move to the next iteration
        y = y0;
        if strcmpi(kernel.type, 'exp2')
            [kernel, ~] = update_kernel_exp2(y(1:T), s, kernel);
        else
            break;
        end
        f0 = f1;
        iter = iter + 1;
        %             return;
    end
end

if debug_on
    avi_file.close();
    close(gcf);
end
end

%% estimate the noise power
function [sn,psdx,ff] = get_noise_fft(Y,options)
% Written by:
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015
% with minor adaption by Pengcheng Zhou, Carnegie Mellon University, 2015
options.noise_range = [.25, .5];
range_ff = options.noise_range;
options.noise_method = 'logmexp';
method = options.noise_method;
options.block_size = [64, 64];
block_size = options.block_size;
options.split_data = false;
split_data = options.split_data;
options.max_timesteps = 3000;

dims = ndims(Y);
sizY = size(Y);
N = min(sizY(end),options.max_timesteps);
if N < sizY(end)
    %Y = reshape(Y,prod(sizY(1:end-1)),[]);
    Y(prod(sizY(1:end-1))*N+1:end) = [];
    Y = reshape(Y,[sizY(1:end-1),N]);
end

Fs = 1;
ff = 0:Fs/N:Fs/2;
indf=ff>range_ff(1);
indf(ff>range_ff(2))=0;
if dims > 1
    d = prod(sizY(1:dims-1));
    Y = reshape(Y,d,N);
    Nb = prod(block_size);
    SN = cell(ceil(d/Nb),1);
    PSDX = cell(ceil(d/Nb),1);
    if ~split_data
        for ind = 1:ceil(d/Nb);
            xdft = fft(Y((ind-1)*Nb+1:min(ind*Nb,d),:),[],2);
            xdft = xdft(:,1: floor(N/2)+1); % FN: floor added.
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(:,2:end-1) = 2*psdx(:,2:end-1);
            %SN{ind} = mean_psd(psdx(:,indf),method);
            switch method
                case 'mean'
                    SN{ind}=sqrt(mean(psdx(:,indf)/2,2));
                case 'median'
                    SN{ind}=sqrt(median(psdx(:,indf)/2),2);
                case 'logmexp'
                    SN{ind} = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
            end
            PSDX{ind} = psdx;
        end
    else
        nc = ceil(d/Nb);
        Yc = mat2cell(Y,[Nb*ones(nc-1,1);d-(nc-1)*Nb],N);
        parfor ind = 1:ceil(d/Nb);
            xdft = fft(Yc{ind},[],2);
            xdft = xdft(:,1:floor(N/2)+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(:,2:end-1) = 2*psdx(:,2:end-1);
            Yc{ind} = [];
            switch method
                case 'mean'
                    SN{ind}=sqrt(mean(psdx(:,indf)/2,2));
                case 'median'
                    SN{ind}=sqrt(median(psdx(:,indf)/2),2);
                case 'logmexp'
                    SN{ind} = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
            end
            
        end
    end
    sn = cell2mat(SN);
else
    xdft = fft(Y);
    xdft = xdft(:,1:floor(N/2)+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(:,2:end-1) = 2*psdx(:,2:end-1);
    switch method
        case 'mean'
            sn = sqrt(mean(psdx(:,indf)/2,2));
        case 'median'
            sn = sqrt(median(psdx(:,indf)/2),2);
        case 'logmexp'
            sn = sqrt(exp(mean(log(psdx(:,indf)/2),2)));
    end
end
psdx = cell2mat(PSDX);
if dims > 2
    sn = reshape(sn,sizY(1:dims-1));
end
end
