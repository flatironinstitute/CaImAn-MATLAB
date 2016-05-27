function [Ain, Cin, bin, fin, center] = initialize_components(Y, K, tau, options)

% Initalize components using a greedy approach followed by hierarchical
% alternative least squares (HALS) NMF. Optional use of spatio-temporal
% downsampling to boost speed.

%Input:
%Y          d1 x d2 x T movie, raw data
%K          number of neurons to extract (default value: 30)
%tau        standard deviation of neuron size (default value: 5)

%options    fine-tuning parameters (optional)
%           options.init_method: method of initialization ('greedy','sparse_NMF','both')
%           options.nIter: number of iterations for shape tuning (default 5)
%           options.gSiz: size of kernel (default 2*tau + 1)
%           options.ssub: spatial downsampling factor (default 1)
%           options.tsub: temporal downsampling factor (default 1)
%           options.nb: rank of background component (default 1)
%           options.save_memory: flag for processing data in chunks to save memory (default 0)
%           options.windowSiz: size of spatial window when computing the median (default 32 x 32)
%           options.chunkSiz: number of timesteps to be processed simultaneously if on save_memory mode (default: 100)
%           options.med_app: number of timesteps to be interleaved for fast (approximate) median calculation (default: 1, no approximation)

%
%Output:
%Ain        (d1*d2) x K matrix, location of each neuron
%Cin        T x K matrix, calcium activity of each neuron
%center     K x 2 matrix, inferred center of each neuron
%bin        (d1*d2) X nb matrix, initialization of spatial background
%fin        nb X T matrix, initalization of temporal background
%res        d1 x d2 x T movie, residual
%
%Authors: Eftychios A. Pnevmatikakis and Pengchen Zhou

if nargin < 4 || isempty(options); options = CNMFSetParms; end
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
% downsample the data

if ~isfield(options, 'ssub'); options.ssub = 1; end; ssub = options.ssub;
if ssub == 1; fprintf('No spatial downsampling is performed. Consider spatial downsampling if the field of view is very large. \n'); end
if ~isfield(options, 'tsub'), options.tsub = 1; end; tsub = options.tsub;
if tsub == 1; fprintf('No temporal downsampling is performed. Consider temporal downsampling if the recording is very long. \n'); end

ndimsY = ndims(Y)-1;
sY = size(Y);
d = sY(1:ndimsY);
T = sY(end);

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
    fprintf('Initializing components with greedy method \n');
    [Ain, Cin, bin, fin] = greedyROI(Y_ds, K, options);
elseif strcmpi(options.init_method, 'greedy_corr')
    fprintf('Initializing components with greedy_corr method \n');
    [Ain, Cin, bin, fin] = greedyROI_corr(Y_ds, K, options);
elseif strcmpi(options.init_method,'sparse_NMF')
    % run sparse_NMF method
    fprintf('Initializing components with sparse NMF \n');
    [Ain,Cin,bin,fin] = sparse_NMF_initialization(Y_ds,K,options_ds);
else
    error('Unknown initialization method')
end

% refine with HALS
fprintf('Refining the initial estimations with HALS...');
[Ain, Cin, bin, fin] = HALS(Y_ds, full(Ain), Cin, bin, fin, options);
K = size(Cin, 1); 
fprintf('  done \n');
%% upsample Ain, Cin, bin, fin
if ndimsY == 2; center = ssub*com(Ain,ds(1),ds(2)); else center = ssub*com(Ain,ds(1),ds(2),ds(3)); end
%Ain = imresize(reshape(Ain, ds(1), ds(2), sum(K)), d);
%Ain = imresize(reshape(full(Ain), [ds, sum(K)]), d);
Ain = imresize(reshape(full(Ain), [ds(1),ds(2), sum(K)*prod(ds)/ds(1)/ds(2)]),[d(1),d(2)]); %,prod(d)/d(1)/d(2)*sum(K)]);
Ain = sparse(reshape(Ain, prod(d), []));
%bin_temp = reshape(bin, ds(1), ds(2), options.nb);
%bin = zeros(d(1),d(2),options.nb);
bin = imresize(reshape(bin,[ds(1),ds(2), options.nb*prod(ds)/ds(1)/ds(2)]),[d(1),d(2)]);
bin = reshape(bin,prod(d),[]);
% bin_temp = reshape(bin, [ds, options.nb]);
% bin = zeros([d,options.nb]);
% for i = 1:options.nb
%     bin(:,:,i) = imresize(bin_temp(:,:,i), d);
% end
%bin = reshape(bin, [], options.nb);
Cin = imresize(Cin, [sum(K), Ts*tsub]);
fin = imresize(fin, [options.nb, Ts*tsub]);
if T ~= Ts*tsub
    Cin = padarray(Cin, [0, T-Ts*tsub], 'post');
    fin = padarray(fin, [0, T-Ts*tsub], fin(end), 'post');
end