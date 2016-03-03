function [P,Y] = preprocess_data(Y,p,options)

% data pre-processing for:
% (i)   identifying and interpolating missing entries (assumed to have the
%       value NaN). Interpolated entries are passed back to Y.
% (ii)  identifying saturated pixels
% (iii) estimating noise level for every pixel
% (iv)  estimating global discrete time constants (if needed)
% This function replaces arpfit, present in the previous versions of the code.

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2015

defoptions.noise_range = [0.25,0.5];            % frequency range over which to estimate the noise
defoptions.noise_method = 'logmexp';            % method for which to estimate the noise level
defoptions.block_size = [64,64];
defoptions.flag_g = false;                          % compute global AR coefficients
defoptions.lags = 5;                                 % number of extra lags when computing the AR coefficients
defoptions.include_noise = 0;                        % include early lags when computing AR coefs
defoptions.split_data = 0;                         % split data into patches for memory reasons

if nargin < 3 || isempty(options); options = defoptions; end
if nargin < 2 || isempty(p); p = 2; end
P.p = p;

if ~isfield(options,'noise_range'); options.noise_range = defoptions.noise_range; end
if ~isfield(options,'noise_method'); options.noise_method = defoptions.noise_method; end
if ~isfield(options,'block_size'); options.block_size = defoptions.block_size; end
if ~isfield(options,'flag_g'); options.flag_g = defoptions.flag_g; end
if ~isfield(options,'lags'); options.lags = defoptions.lags; end
if ~isfield(options,'include_noise'); options.include_noise = defoptions.include_noise; end; include_noise = options.include_noise;
if ~isfield(options,'split_data'); split_data = defoptions.split_data; else split_data = options.split_data; end

%% interpolate missing data

if any(isnan(Y(:)))
    Y_interp = interp_missing_data(Y);      % interpolate missing data
    mis_data = find(Y_interp);
    Y(mis_data) = Y_interp(mis_data);       % introduce interpolated values for initialization
else
    Y_interp = sparse(size(Y));
    mis_data = [];
end
P.mis_values = full(Y_interp(mis_data));
P.mis_entries = mis_data;

%% indentify saturated pixels

P.pixels = find_unsaturatedPixels(Y);                % pixels that do not exhibit saturation

%% estimate noise levels

fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
[sn,psx] = get_noise_fft(Y,options);
P.sn = sn(:);
fprintf('  done \n');

%% cluster pixels based on PSD

psdx = sqrt(psx(:,3:end));
X = psdx(:,1:min(size(psdx,2),1500));
X = bsxfun(@minus,X,mean(X,2));     % center
X = spdiags(std(X,[],2),0,size(X,1),size(X,1))\X;
[L,Cx] = kmeans_pp(X',2);
[~,ind] = min(sum(Cx(end-49:end,:),1));
P.active_pixels = (L==ind);

% [P.W,P.H] = nnmf(sqrt(psdx(:,3:end)),2); %,'h0',H0);
% r = sort(rand(1,size(psdx,2)-2),'descend');
% H = [r/norm(r); ones(1,length(r))/sqrt(length(r))];
% for iter = 1:100
%     W = max((H*H')\(H*psdx'),0)';
%     H = max((W'*W)\(W'*psdx),0);
% end
% P.W = W;
% P.H = H;

%% estimate global time constants

if options.flag_g
    if ndims(Y) == 3
        Y = reshape(Y,size(Y,1)*size(Y,2),size(Y,3));
    end

    ff = options.pixels;
    np = length(ff);

    fprintf('Estimating time constant through autocorrelation function.. \n');
    tt1 = tic;
    mp = max(p);
    lags = options.lags + mp;
    if split_data 
        Ycl = mat2cell(Y(ff,:),ones(np,1),size(Y,2));
        XC = cell(np,1);
        parfor j = 1:np
            XC{j} = xcov(Ycl{j},lags,'biased');
        end
        XC = cell2mat(XC);   
    else
        XC = zeros(np,2*lags+1);
        for j = 1:np
            XC(j,:) = xcov(Y(ff(j),:),lags,'biased');
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
    P.g = ph(:);
end