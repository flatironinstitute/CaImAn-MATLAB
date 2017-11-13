function [fitness,erfc,sd_r,md] = compute_event_exceptionality(traces,N,robust_std)
%{
    Define a metric and order components according to the probabilty if some "exceptional events" (like a spike).
    Suvh probability is defined as the likeihood of observing the actual trace value over N samples given an estimated noise distribution.
    The function first estimates the noise distribution by considering the dispersion around the mode.
    This is done only using values lower than the mode. The estimation of the noise std is made robust by using the approximation std=iqr/1.349.
    Then, the probavility of having N consecutive eventsis estimated.
    This probability is used to order the components.
 
    
    Parameters:
    -----------

    traces: array
        Fluorescence traces

 
    Returns:
    --------

    fitness: array
        value estimate of the quality of components (the lesser the better)

    erfc: array
        probability at each time step of observing the N consequtive actual trace values given the distribution of noise

    noise_est: array
        the components ordered according to the fitness

%}

T=size(traces,2);

md = max(mode_robust(traces, 2),0);

if ~exist('N','var'); N=5; end
if ~exist('robust_std','var')
    robust_std = false;
end

ff1 = bsxfun(@minus,traces,md');
% only consider values under the mode to determine the noise standard deviation
ff1 = -ff1 .* (ff1 < 0);
if robust_std
    % compute 25 percentile
    ff1 = sort(ff1, 2);
    ff1(ff1 == 0) = nan;
    Ns = round(sum(ff1 > 0, 2) * .5);
    iqr_h = zeros(size(traces,1));
    idx = 1;
    for idx = 1:size(ff1,1) 
        el = ff1(idx,:);
        iqr_h(idx) = ff1(idx, -Ns(idx));
    end
    % approximate standard deviation as iqr/1.349
    sd_r = 2 * iqr_h / 1.349;
    
else
    Ns = sum(ff1 > 0, 2);
    sd_r = sqrt(sum(ff1.^2, 2)./ Ns);
end
% compute z value
z = bsxfun(@times,bsxfun(@minus,traces,md'),1./(3 * sd_r));
% probability of observing values larger or equal to z given normal
% distribution with mean md and std sd_r
mu = 0;
sigma = 1;
pd = makedist('Normal',mu,sigma);
erf = 1 - cdf(pd,z);
% use logarithm so that multiplication becomes sum
erf = log(erf);
filt = ones(1,N);
% moving sum
erfc = conv2(1,filt,erf,'same');
erfc = erfc(:,1:T);

% select the maximum value of such probability for each trace
fitness = min(erfc, [], 2);

% ordered = np.argsort(fitness)





function outp = hsm_(data)
outp = [];
if numel(data) == 1
    outp = data(1);
elseif numel(data) == 2
    outp = mean(data(:));
elseif numel(data) == 3
    i1 = data(2) - data(1);
    i2 = data(3) - data(2);
    if i1 < i2
        outp = mean(data(1:2));
    elseif i2 > i1
        outp = mean(data(2:end));
    else
        outp = data(2);
    end
else
    
    % wMin = data[-1] - data[0]
    wMin = inf;
    N = idivide(int32(numel(data)),2,'floor') + mod(numel(data),2)-1;
    
    for i = 1:N
        w = data(i + N - 1) - data(i);
        if w < wMin
            wMin = w;
            j = i;
        end
    end
    outp = hsm_(data(j:j + N));
end
                
function dataMode = mode_robust(inputData, axis)
    %{
    Robust estimator of the mode of a data set using the half-sample mode.
    
    .. versionadded: 1.0.3
    """
    %}
    if ~isempty(axis)
        
        if axis == 2
            dataMode = arrayfun(@(n)  mode_robust(inputData(n,:),[]), 1:size(inputData,1));
        elseif axis == 1
            dataMode = arrayfun(@(n)  mode_robust(inputData(:,n),[]), 1:size(inputData,2));
        else
            error('axis can be 1 or two')
        end
    else
        % Create the function that we can use for the half-sample mode
        data = inputData(:);
        % The data need to be sorted for this to work
        data = sort(data);
        % Find the mode
        dataMode = hsm_(data);
    end
    