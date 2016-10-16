function [Ain, Cin, bin, fin, center] = initialize_components(Y, K, tau, options, P)

% Initalize components using a greedy approach followed by hierarchical
% alternative least squares (HALS) NMF. Optional use of spatio-temporal
% downsampling to boost speed.

%Input:
%Y          d1 x d2 x T movie, raw data
%K          number of neurons to extract (default value: 30)
%tau        standard deviation of neuron size (default value: 5)

%options    fine-tuning parameters (optional)
%           options.init_method: method of initialization ('greedy','sparse_NMF','HALS')
%           options.
%           options.nIter: number of iterations for shape tuning (default 5)
%           options.gSiz: size of kernel (default 2*tau + 1)
%           options.ssub: spatial downsampling factor (default 1)
%           options.tsub: temporal downsampling factor (default 1)
%           options.nb: rank of background component (default 1)
%           options.save_memory: flag for processing data in chunks to save memory (default 0)
%           options.windowSiz: size of spatial window when computing the median (default 32 x 32)
%           options.chunkSiz: number of timesteps to be processed simultaneously if on save_memory mode (default: 100)
%           options.med_app: number of timesteps to be interleaved for fast (approximate) median calculation (default: 1, no approximation)
%           options.rem_prct: percentile to be removed before initialization (default: 20)

% P         parameter struct used for normalization by noise and user feed component centroids (optional)

%
%Output:
%Ain        (d1*d2) x K matrix, location of each neuron
%Cin        T x K matrix, calcium activity of each neuron
%center     K x 2 matrix, inferred center of each neuron
%bin        (d1*d2) X nb matrix, initialization of spatial background
%fin        nb X T matrix, initalization of temporal background
%res        d1 x d2 x T movie, residual
%
%Authors: Eftychios A. Pnevmatikakis and Pengchen Zhou, with inputs from Weijian Yang


defoptions = CNMFSetParms;
if nargin < 4 || isempty(options); options = defoptions; end
if nargin < 2 || isempty(K)
    K = 30;
    fprintf('Number of components to be detected not specified. Using the default value 30. \n');
end
if nargin < 3 || isempty(tau)
    options.gSig = 5;
    fprintf('Standard deviation for neuron not specified. Using the default value 5. \n');
else options.gSig = tau;
end

if ~isfield(options,'init_method'); options.init_method = 'greedy'; end
if ~isfield(options,'rem_prct') || isempty(options.rem_prct); options.rem_prct = defoptions.rem_prct; end
% downsample the data

if ~isfield(options, 'ssub'); options.ssub = 1; end; ssub = options.ssub;
if ssub == 1; fprintf('No spatial downsampling is performed. Consider spatial downsampling if the field of view is very large. \n'); end
if ~isfield(options, 'tsub'), options.tsub = 1; end; tsub = options.tsub;
if tsub == 1; fprintf('No temporal downsampling is performed. Consider temporal downsampling if the recording is very long. \n'); end

if ~isfield(options,'noise_norm') || isempty(options.noise_norm)
    options.noise_norm = defoptions.noise_norm; % normalization by noise (true if P is present)
end

if nargin < 5
    if options.noise_norm
        warning('Normalization by noise value is not performed since noise values are not provided. \n');
    end
    options.noise_norm = false;
end

ndimsY = ndims(Y)-1;
sY = size(Y);
d = sY(1:ndimsY);
T = sY(end);

if options.noise_norm
    mY = mean(Y,ndims(Y));
    norm_image = mY + median(mY(:)) + 1e-4;
    min_noise = norm_image;
    %min_noise = prctile(P.sn(P.sn>0),options.noise_norm_prctile);
    %Y = bsxfun(@times,Y,reshape(1./max(P.sn,min_noise),d));
    Y = bsxfun(@times,Y,reshape(1./double(min_noise),d));
end

ds = d;
ds(1:2) = ceil(d(1:2)/ssub); % do not subsample along z axis
%d1s = ceil(d1/ssub);        %size of downsampled image
%d2s = ceil(d2/ssub);
Ts = floor(T/tsub);         %reduced number of frames
% spatial downsampling
fprintf('starting resampling \n')
if ssub~=1;
    if ndimsY == 2; Y_ds = imresize(Y, [ds(1), ds(2)], 'box'); end
    if ndimsY == 3;
        Y_ds = zeros([ds(1:2),T,ds(end)]);
        for z = 1:ds(3)
            Y_ds(:,:,:,z) = imresize(squeeze(Y(:,:,z,:)), [ds(1), ds(2)], 'box');
        end
        Y_ds = permute(Y_ds,[1,2,4,3]);
    end
else
    Y_ds = Y;
end
% temporal downsampling
if tsub~=1
    if ndimsY == 2; Y_ds = squeeze(mean(reshape(Y_ds(:, :, 1:(Ts*tsub)),ds(1), ds(2), tsub, Ts), 3)); end
    if ndimsY == 3; Y_ds = squeeze(mean(reshape(Y_ds(:, :, :, 1:(Ts*tsub)),ds(1), ds(2), ds(3), tsub, Ts), 4)); end
end

options_ds = options;
options_ds.d1 = ds(1);
options_ds.d2 = ds(2);

if strcmpi(options.init_method,'greedy')
    % run greedy method
    if nargin < 5 || ~isfield(P,'ROI_list')
        ROI_list = [];
    else
        ROI_list = round(P.ROI_list/ssub);
        K = size(ROI_list,1);
    end
    fprintf('Initializing components with greedy method \n');
    [Ain, Cin, bin, fin] = greedyROI(Y_ds, K, options, ROI_list);
elseif strcmpi(options.init_method, 'greedy_corr')
    fprintf('Initializing components with greedy_corr method \n');
    [Ain, Cin, bin, fin] = greedyROI_corr(Y_ds, K, options);
elseif strcmpi(options.init_method,'sparse_NMF')
    % run sparse_NMF method
    fprintf('Initializing components with sparse NMF \n');
    [Ain,Cin,bin,fin] = sparse_NMF_initialization(Y_ds,K,options_ds);
elseif strcmpi(options.init_method,'HALS')
    fprintf('Initializing components with HALS \n');
    [Ain,Cin,bin,fin] = HALS_initialization(Y_ds,K,options_ds);
else
    error('Unknown initialization method')
end

% refine with HALS
fprintf('Refining initial estimates with HALS...');
[Ain, Cin, bin, fin] = HALS(Y_ds, full(Ain), Cin, bin, fin, options_ds);
fprintf('  done \n');
%% upsample Ain, Cin, bin, fin
if nargout == 5
    if ndimsY == 2; center = ssub*com(Ain,ds(1),ds(2)); else center = ssub*com(Ain,ds(1),ds(2),ds(3)); end
end

Ain = imresize(reshape(full(Ain), [ds(1),ds(2), size(Ain,2)*prod(ds)/ds(1)/ds(2)]),[d(1),d(2)]); %,prod(d)/d(1)/d(2)*sum(K)]);
Ain = max(sparse(reshape(Ain, prod(d), [])),0);

bin = imresize(reshape(bin,[ds(1),ds(2), options.nb*prod(ds)/ds(1)/ds(2)]),[d(1),d(2)]);
bin = max(double(reshape(bin,prod(d),[])),0);

if options.noise_norm
    %Ain = bsxfun(@times,Ain,double(max(P.sn(:),min_noise)));
    %bin = bsxfun(@times,bin,max(P.sn(:),min_noise));
    Ain = bsxfun(@times,Ain,double(min_noise(:)));
    bin = bsxfun(@times,bin,double(min_noise(:)));
end
Cin = max(imresize(Cin, [size(Cin, 1), Ts*tsub]),0);
fin = max(imresize(fin, [options.nb, Ts*tsub]),0);
if T ~= Ts*tsub
    Cin = padarray(Cin, [0, T-Ts*tsub], 'post');
    fin = padarray(fin, [0, T-Ts*tsub], fin(end), 'post');
end
