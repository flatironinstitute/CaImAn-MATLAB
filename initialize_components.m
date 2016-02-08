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

[d1,d2,T] = size(Y);
d1s = ceil(d1/ssub);        %size of downsampled image
d2s = ceil(d2/ssub);
Ts = floor(T/tsub);         %reduced number of frames
% spatial downsampling
if ssub~=1, Y_ds = imresize(Y, [d1s, d2s], 'box'); else Y_ds = Y; end
% temporal downsampling
if tsub~=1
    Y_ds = squeeze(mean(reshape(Y_ds(:, :, 1:(Ts*tsub)),d1s, d2s, tsub, Ts), 3));
end

if strcmpi(options.init_method,'greedy')
    % run greedy method
    fprintf('Initializing components with greedy method \n');
    [Ain, Cin, bin, fin] = greedyROI2d(Y_ds, K, options);
elseif strcmpi(options.init_method,'sparse_NMF')
    % run sparse_NMF method
    fprintf('Initializing components with sparse NMF \n');
    [Ain,Cin,bin,fin] = sparse_NMF_initialization(Y_ds,K,options);
else
    error('Unknown initialization method')
end

% refine with HALS
fprintf('Refining initial estimates with HALS \n');
[Ain, Cin, bin, fin] = HALS_2d(Y_ds, full(Ain), Cin, bin, fin, options); 

%% upsample Ain, Cin, bin, fin
center = ssub*com(Ain,d1s,d2s); 
Ain = imresize(reshape(Ain, d1s, d2s, sum(K)), [d1, d2]);
Ain = reshape(Ain, d1*d2, []); 
bin_temp = reshape(bin, d1s, d2s, options.nb);
bin = zeros(d1,d2,options.nb);
for i = 1:options.nb
    bin(:,:,i) = imresize(bin_temp(:,:,i), [d1, d2]);
end
bin = reshape(bin, [], options.nb); 
Cin = imresize(Cin, [sum(K), Ts*tsub]);
fin = imresize(fin, [options.nb, Ts*tsub]);
if T ~= Ts*tsub
    Cin = padarray(Cin, [0, T-Ts*tsub], 'post'); 
    fin = padarray(fin, [0, T-Ts*tsub], fin(end), 'post');   
end