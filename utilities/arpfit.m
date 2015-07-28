function P = arpfit(Y,p,options,sn)

% Estimate a global time constant of order p and noise power for each pixel
% The time constants are estimated by the decay of the autorcorrelation
% function. The noise is estimated from the tails of the power spectral density.

defoptions.noise_range = [0.25,0.5];            % frequency range over which to estimate the noise
defoptions.noise_method = 'logmexp';            % method for which to estimate the noise level
defoptions.block_size = [64,64];
defoptions.lags = 5;                                 % number of extra lags when computing the AR coefficients
defoptions.include_noise = 0;                        % include early lags when computing AR coefs
defoptions.pixels = 1:numel(Y)/size(Y,ndims(Y));     % pixels to include when computing the AR coefs
defoptions.use_parallel = 0;                         % split data into patches for memory reasons

if nargin < 4
    sn = [];
    if nargin < 3 || isempty(options)
        options = defoptions;
    end
end

if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'include_noise'); options.include_noise = defoptions.include_noise; end; include_noise = options.include_noise;
if ~isfield(options,'pixels'); options.pixels = defoptions.pixels; end
if ~isfield(options,'use_parallel'); use_parallel = defoptions.use_parallel; else use_parallel = options.use_parallel; end

if isempty(sn)
    fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
    sn = get_noise_fft(Y,options.noise_range,options.noise_method,options.block_size);
    fprintf('  done \n');
end

find_bas = 0;   % estimate baseline from quantiles

if ndims(Y) == 3
    Y = reshape(Y,size(Y,1)*size(Y,2),size(Y,3));
end

ff = options.pixels;
np = length(ff);

fprintf('Estimating time constant through autocorrelation function.. \n');
tt1 = tic;
mp = max(p);
lags = options.lags + mp;
if use_parallel 
    Ycl = mat2cell(Y(ff,:),ones(np,1),size(Y,2));
    XC = cell(np,1);
    parfor j = 1:np
        XC{j} = xcov(Ycl{j},lags,'biased');
    end
    XC = cell2mat(XC);   
else
    XC = zeros(np,2*lags+1);
    for j = 1:np
        XC(j,:) = xcov(Y(j,:),lags,'biased');
    end
end
    
    gv = zeros(np*lags,1);
    if ~include_noise
        g =  XC(:,lags:-1:1);
        clear XC;
        lags = lags - p;
    end
    A = zeros(np*lags,p);
    for i = 1:np
        if ~include_noise
            A((i-1)*lags + (1:lags),:) = toeplitz(g(i,p:p+lags-1),g(i,p:-1:1));
        else
            A((i-1)*lags + (1:lags),:) = toeplitz(XC(i,lags+(1:lags)),XC(i,lags+(1:p))) - sn(i)^2*eye(lags,p);
            gv((i-1)*lags + (1:lags)) = XC(i,lags+2:end)';
        end
    end
    if ~include_noise
        gv = g(:,p+1:end)';
    end
    %ph = pinv(A)*gv(:);
    ph = A\gv(:);
    disp(ph);
    fprintf('Done after %2.2f seconds. \n',toc(tt1));


P.sn = sn(:);
P.g = ph(:);


    function [sn,psdx,ff] = get_noise_fft(Y,range_ff,method,block_size)
        
        dims = ndims(Y);
        sizY = size(Y);
        N = sizY(end);
        Fs = 1;        
        ff = 0:Fs/N:Fs/2;
        indf=ff>range_ff(1);
        indf(ff>range_ff(2))=0;
        if dims > 1
            d = prod(sizY(1:dims-1));
            Y = reshape(Y,d,N);
            Nb = prod(block_size);
            SN = cell(ceil(d/Nb),1);
            if isempty(which('parpool'));
                for ind = 1:ceil(d/Nb); 
                    xdft = fft(Y((ind-1)*Nb+1:min(ind*Nb,d),:),[],2); 
                    xdft = xdft(:,1:N/2+1);
                    psdx = (1/(Fs*N)) * abs(xdft).^2;
                    psdx(:,2:end-1) = 2*psdx(:,2:end-1);
                    SN{ind} = mean_psd(psdx(:,indf),method);
                end
            else
                nc = ceil(d/Nb);
                Yc = mat2cell(Y,[Nb*ones(nc-1,1);d-(nc-1)*Nb],N);
                parfor ind = 1:ceil(d/Nb); 
                    xdft = fft(Yc{ind},[],2); 
                    xdft = xdft(:,1:N/2+1);
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
            xdft = xdft(:,1:N/2+1);
            psdx = (1/(Fs*N)) * abs(xdft).^2;
            psdx(:,2:end-1) = 2*psdx(:,2:end-1);
            sn = mean_psd(psdx(:,indf),'method');
        end
        clear psdx
        if dims > 2
            sn = reshape(sn,sizY(1:dims-1));
        end

    end

    function mp = mean_psd(y,method)
        switch method
            case 'mean'
                mp=sqrt(mean(y/2,2));
            case 'median'
                mp=sqrt(median(y/2),2);
            case 'logmexp'
                mp = sqrt(exp(mean(log(y/2),2)));
        end
    end

end