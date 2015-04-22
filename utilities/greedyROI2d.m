function [basis, trace, center, data] = greedyROI2d(data, K, params)
%greedyROI2d using greedy algorithm to identify neurons in 2d calcium movie
%
%Usage:     [basis, trace, center, res] = greedyROI2d(data, K, params)
%
%Input:
%data       M x N x T movie, raw data, each column is a vectorized movie
%K          number of neurons to extract
%params     tuning parameter for fine-tuning the shape (optional)
%           params.nIter: number of iterations for shape tuning (default 5)
%           params.gSig: variance of Gaussian kernel to use (default 5)
%           params.gSiz: size of kernel (default 41)
%
%Output:
%basis      M x N x K matrix, location of each neuron
%trace      T x K matrix, calcium activity of each neuron
%center     K x 2 matrix, inferred center of each neuron
%res        M x N x T movie, residual
%
%Authur: Yuanjun Gao

[M, N, T] = size(data);

med = median(data, 3);
data = bsxfun(@minus, data, med);

if ~exist('K', 'var'),  %find K neurons
    K = 30;  
    warning(['number of neurons are not specified, set to be the default value', num2str(K)]);
end

if ~exist('params', 'var') params = []; end

if ~isfield(params, 'gSig'), gSig = [5, 5]; 
elseif length(params.gSig) == 1, gSig = params.gSig + zeros(1,2); end

if ~isfield(params, 'gSiz'), gSiz = [11, 11];
elseif length(params.gSiz) == 1, gSiz = params.gSiz + zeros(1,2); end

if ~isfield(params, 'nIter'), nIter = 5; 
else nIter = params.nIter; end

basis = zeros(M, N, K);
trace = zeros(T, K);
center = zeros(K, 2);

gHalf = floor(gSiz / 2); %half size of the kernel, used to calculate margin
gSiz = 2 * gHalf + 1; %actual size

%scan the whole image (only need to do this at the first iteration)
rho = imblur(data, gSig, gSiz); %covariance of data and basis
v = sum(rho.^2, 3); %variance explained

for k = 1:K,    
    [~, ind] = max(v(:));
    [iHat, jHat] = ind2sub([M, N], ind);
    center(k, 1) = iHat; center(k, 2) = jHat;
        
    iSig = [max(iHat - gHalf(1), 1), min(iHat + gHalf(1), M)]; iSigLen = iSig(2) - iSig(1) + 1;
    jSig = [max(jHat - gHalf(2), 1), min(jHat + gHalf(2), N)]; jSigLen = jSig(2) - jSig(1) + 1;
        
    %fine tune the shape
    dataTemp = data(iSig(1):iSig(2), jSig(1):jSig(2), :);
    traceTemp = rho(iHat, jHat, :);
    [coef, score] = finetune2d(dataTemp, traceTemp, nIter);        
    
    dataSig = bsxfun(@times, coef, reshape(score, [1,1,T]));
    basis(iSig(1):iSig(2), jSig(1):jSig(2), k) = coef;
    trace(:, k) = score';
            
    data(iSig(1):iSig(2), jSig(1):jSig(2), :) = data(iSig(1):iSig(2), jSig(1):jSig(2), :) - dataSig; %update residual
    fprintf('found %i out of %i neurons..\n', k,K);
    
    %get next basis;
    if k < K,
        iMod = [max(iHat - 2 * gHalf(1), 1), min(iHat + 2 * gHalf(1), M)]; iModLen = iMod(2) - iMod(1) + 1;%patches to modify
        jMod = [max(jHat - 2 * gHalf(2), 1), min(jHat + 2 * gHalf(2), N)]; jModLen = jMod(2) - jMod(1) + 1;
        iLag = iSig - iMod(1) + 1; %relative location of iSig in the small patch
        jLag = jSig - jMod(1) + 1;
        dataTemp = zeros(iModLen, jModLen);
        dataTemp(iLag(1):iLag(2), jLag(1):jLag(2)) = reshape(coef, [iSigLen, jSigLen]);
        dataTemp = imblur(dataTemp, gSig, gSiz, 2);
        rhoTemp = bsxfun(@times, dataTemp, reshape(score, [1,1,T]));
        rhoTemp = rho(iMod(1):iMod(2), jMod(1):jMod(2), :) - rhoTemp;
        rho(iMod(1):iMod(2), jMod(1):jMod(2), :) = rhoTemp;
        v(iMod(1):iMod(2), jMod(1):jMod(2)) = sum(rhoTemp.^2, 3);
    end
end

%basis = reshape(basis, [M * N, K]);


function [basis, trace] = finetune2d(data, trace, nIter)
%using matrix factorization with lasso penalty to fine-tune the basis
%
%Input:
%data   M x N x T matrix, small patch containing one neuron
%trace  initial value for trace
%nIter  number of coordinate descent steps
%
%Output:
%basis  M x N matrix, result of the fine-tuned neuron shape
%trace  1 x T matrix, result of the neuron

T = size(data, 3);
if ~exist('nIter', 'var'), nIter = 1; end

%do block coordinate descent
for iter = 1:nIter,
    %update basis
    a = sum(trace.^2); %scale by num. of observation
    b = sum(bsxfun(@times, data, reshape(trace, [1,1,T])), 3);
    basis = max(b / a, 0);
    basisNorm = norm(basis(:));
    if basisNorm > 0, 
        basis = basis / basisNorm; 
    else
        fprintf('get degenerate basis!\n')
        break
    end
    
    %updating trace
    trace = bsxfun(@times, data, basis);
    trace = squeeze(sum(sum(trace, 1), 2));            
end
end



function data = imblur(data, sig, siz, nDimBlur)
%Gaussian blur for high dimensional data
%Input:
%data       original data
%sig        std of gaussian kernel
%siz        size of kernel
%nDimBlur   number of dims to blur (default: ndims(data) - 1)
%
%Output:
%data       result after the Gaussian blur

if ~exist('nDimBlur', 'var'), nDimBlur = ndims(data) - 1; 
else nDimBlur = min(nDimBlur, ndims(data)); end

if length(sig) == 1, sig = sig * ones(1,nDimBlur); end
assert(nDimBlur == length(sig));

if length(siz) == 1, siz = siz * ones(1,nDimBlur); end
assert(nDimBlur == length(siz));

for i = 1:nDimBlur,
    if sig(i) > 0,
        x = -floor(siz(i) / 2):floor(siz(i) / 2);
        H = exp(-(x.^2/ (2 * sig(i)^2)));
        H = H' / norm(H(:));
        if nDimBlur > 1,
            indH = 1:nDimBlur; indH(i) = 1; indH(1) = i;
            H = permute(H, indH);
        end
        data = imfilter(data, H, 'same', 0);
    end
end
end


end


