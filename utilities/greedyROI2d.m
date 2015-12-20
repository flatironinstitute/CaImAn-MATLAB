function [Ain, Cin, bin, fin, center, res] = greedyROI2d(data, K, params)
%greedyROI2d using greedy algorithm to identify neurons in 2d calcium movie
%
%Usage:     [basis, trace, center, res] = greedyROI2d(data, K, params)
%
%Input:
%data       M x N x T movie, raw data, each column is a vectorized movie
%K          number of neurons to extract (if K is a vector then the algorithm is run multiple times with different parameters)
%params     tuning parameter for fine-tuning the shape (optional)
%           params.nIter: number of iterations for shape tuning (default 5)
%           params.gSig: variance of Gaussian kernel to use (default 5) If params.gSig is a cell, then the algorithm is run with multiple times with different parameters
%           params.gSiz: size of kernel (default 2*gSiz+1)
%           params.nb: rank of background component (default 1)
%           params.save_memory: flag for processing data in chunks to save memory (default 0)
%           params.windowSiz: size of spatial window when computing the median (default 32 x 32)
%           params.chunkSiz: number of timesteps to be processed simultaneously if on save_memory mode (default: 100)
%           params.med_app: number of timesteps to be interleaved for fast (approximate) median calculation (default: 1, no approximation)



%Output:
%Ain        (MN) x K matrix, location of each neuron
%Cin        T x K matrix, calcium activity of each neuron
%center     K x 2 matrix, inferred center of each neuron
%bin        (MN) X nb matrix, initialization of spatial background
%fin        nb X T matrix, initalization of temporal background
%res        M x N x T movie, residual
%
%Author: Yuanjun Gao with modifications from Eftychios A. Pnevmatikakis

[M, N, T] = size(data);

if ~exist('K', 'var'),  %find K neurons
    K = 30;  
    warning(['number of neurons are not specified, set to be the default value', num2str(K)]);
end

if ~exist('params', 'var') params = []; end

if ~isfield(params, 'gSig') || isempty(params.gSig); params.gSig = [5, 5]; 
elseif length(params.gSig) == 1, params.gSig = params.gSig + zeros(1,2); 
    %else gSig = params.gSig; 
end

if ~isfield(params, 'gSiz') || isempty(params.gSiz); params.gSiz = 2*params.gSig + 1;
elseif length(params.gSiz) == 1, params.gSiz = params.gSiz + zeros(1,2); 
    %else gSiz = params.gSiz; 
end

if isfield(params,'ssub'); params.gSig = params.gSig/params.ssub; params.gSiz = ceil(params.gSiz/params.ssub); end

if ~isfield(params,'nb'), nb = 1; else nb = params.nb; end

if ~isfield(params, 'nIter'), nIter = 5; 
else nIter = params.nIter; end

if ~isfield(params, 'save_memory'), save_memory = 0;
else save_memory = params.save_memory; end
    
if ~isfield(params, 'chunkSiz'), chunkSiz = 100;
else chunkSiz = params.chunkSiz; end

if ~isfield(params, 'windowSiz'), windowSiz = 32;
else windowSiz = params.windowSiz; end

if ~isfield(params, 'med_app'), med_app = 1;
else med_app = params.med_app; end

Tint = 1:med_app:T;
if save_memory
    med = zeros(M,N);
    for ii = 1:ceil(M/windowSiz)
        intx = (ii-1)*windowSiz+1:min(ii*windowSiz,M);
        for jj = 1:ceil(N/windowSiz)
            inty = (jj-1)*windowSiz+1:min(jj*windowSiz,N);
            med(intx,inty) = median(data(intx,inty,Tint),3);
        end
    end 
else
    med = median(data(:,:,Tint), 3);
end

data = bsxfun(@minus, data, med);
if iscell(params.gSig); params.gSig = cell2mat(params.gSig); params.gSiz = cell2mat(params.gSiz); end
if length(K) > 1  % order size of components to be found in descending order
    [~,ord] = sort(sum(params.gSig,2),'descend');
    K = K(ord);
    params.gSig = params.gSig(ord,:);
    params.gSiz = params.gSiz(ord,:);
end
Ain = zeros(M*N,sum(K));
Cin = zeros(sum(K),T);
center = zeros(sum(K),2);

for r = 1:length(K)
    gSig = params.gSig(r,:);
    gSiz = params.gSiz(r,:);
    
    gHalf = floor(gSiz / 2); %half size of the kernel, used to calculate margin
    gSiz = 2 * gHalf + 1; %actual size

    basis = zeros(M, N, K(r));
    %basis = spalloc(M*N,K,K*prod(gSiz));
    trace = zeros(T, K(r));
    centers = zeros(K(r), 2);


    %scan the whole image (only need to do this at the first iteration)
    rho = imblur(data, gSig, gSiz, ndims(data)-1, save_memory, chunkSiz); %covariance of data and basis
    v = sum(rho.^2, 3); %variance explained

    for k = 1:K(r),    
        [~, ind] = max(v(:));
        [iHat, jHat] = ind2sub([M, N], ind);
        centers(k, 1) = iHat; centers(k, 2) = jHat;

        iSig = [max(iHat - gHalf(1), 1), min(iHat + gHalf(1), M)]; iSigLen = iSig(2) - iSig(1) + 1;
        jSig = [max(jHat - gHalf(2), 1), min(jHat + gHalf(2), N)]; jSigLen = jSig(2) - jSig(1) + 1;

        %fine tune the shape
        dataTemp = data(iSig(1):iSig(2), jSig(1):jSig(2), :);
        traceTemp = rho(iHat, jHat, :);
        [coef, score] = finetune2d(dataTemp, traceTemp, nIter);        

        dataSig = bsxfun(@times, coef, reshape(score, [1,1,T]));
        %[xSig,ySig] = meshgrid(iSig(1):iSig(2),jSig(1):jSig(2));
        basis(iSig(1):iSig(2), jSig(1):jSig(2), k) = coef;
        %basis(sub2ind([M,N],xSig(:),ySig(:)),k) = coef(:);
        trace(:, k) = score';

        data(iSig(1):iSig(2), jSig(1):jSig(2), :) = data(iSig(1):iSig(2), jSig(1):jSig(2), :) - dataSig; %update residual
        if mod(k,10) == 0; fprintf('found %i out of %i neurons..\n', k,K(r)); end

        %get next basis;
        if k < K(r)
            iMod = [max(iHat - 2 * gHalf(1), 1), min(iHat + 2 * gHalf(1), M)]; iModLen = iMod(2) - iMod(1) + 1;%patches to modify
            jMod = [max(jHat - 2 * gHalf(2), 1), min(jHat + 2 * gHalf(2), N)]; jModLen = jMod(2) - jMod(1) + 1;
            iLag = iSig - iMod(1) + 1; %relative location of iSig in the small patch
            jLag = jSig - jMod(1) + 1;
            dataTemp = zeros(iModLen, jModLen);
            dataTemp(iLag(1):iLag(2), jLag(1):jLag(2)) = reshape(coef, [iSigLen, jSigLen]);
            dataTemp = imblur(dataTemp, gSig, gSiz, 2, 0, chunkSiz);
            rhoTemp = bsxfun(@times, dataTemp, reshape(score, [1,1,T]));
            rhoTemp = rho(iMod(1):iMod(2), jMod(1):jMod(2), :) - rhoTemp;
            rho(iMod(1):iMod(2), jMod(1):jMod(2), :) = rhoTemp;
            v(iMod(1):iMod(2), jMod(1):jMod(2)) = sum(rhoTemp.^2, 3);
        end
    end
    Ain(:,sum(K(1:r-1))+1:sum(K(1:r))) = sparse(reshape(basis,M*N,K(r)));  
    Cin(sum(K(1:r-1))+1:sum(K(1:r)),:) = trace';
    center(sum(K(1:r-1))+1:sum(K(1:r)),:) = centers;
end
res = reshape(data,M*N,T) + repmat(med(:),1,T);
[bin,fin] = nnmf(max(res,0),nb);


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

function data = imblur(data, sig, siz, nDimBlur, save_memory, chunkSiz)
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
        if save_memory
            L = size(data,ndims(data));
            for ci = 1:ceil(L/chunkSiz)
                int = (ci-1)*chunkSiz+1:min(ci*chunkSiz,L);
                data(:,:,int) = imfilter(data(:,:,int), H, 'same', 0);
            end
        else
            data = imfilter(data, H, 'same', 0);
        end
    end
end
end

end
