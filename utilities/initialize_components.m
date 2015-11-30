function [Ain, Cin, bin, fin, center] = initialize_components(Y, nr, params)

% Initalize components using a greedy approach followed by hierarchical
% alternative least squares (HALS) NMF. Optional use of spatio-temporal
% downsampling to boost speed.

%Input:
%Y          d1 x d2 x T movie, raw data
%nr         number of neurons to extract
%params     tuning parameter for fine-tuning the shape (optional)
%           params.nIter: number of iterations for shape tuning (default 5)
%           params.gSig: variance of Gaussian kernel to use (default 5)
%           params.gSiz: size of kernel (default 2*params.gSiz + 1)
%           params.ssub: spatial downsampling factor (default 1)
%           params.tsub: temporal downsampling factor (default 1)
%           params.nb: rank of background component (default 1)
%           params.save_memory: flag for processing data in chunks to save memory (default 0)
%           params.windowSiz: size of spatial window when computing the median (default 32 x 32)
%           params.chunkSiz: number of timesteps to be processed simultaneously if on save_memory mode (default: 100)
%           params.med_app: number of timesteps to be interleaved for fast (approximate) median calculation (default: 1, no approximation)

%
%Output:
%Ain        (d1*d2) x nr matrix, location of each neuron
%Cin        T x nr matrix, calcium activity of each neuron
%center     nr x 2 matrix, inferred center of each neuron
%bin        (d1*d2) X nb matrix, initialization of spatial background
%fin        nb X T matrix, initalization of temporal background
%res        d1 x d2 x T movie, residual
%
%Authors: Eftychios A. Pnevmatikakis and Pengchen Zhou

% downsample the data

if isfield(params, 'ssub'), ssub = params.ssub; else ssub = 1; end
if isfield(params, 'tsub'), tsub = params.tsub; else tsub = 1; end
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

% run greedy method
fprintf('Initializing components with greedy method \n');
[Ain, Cin, bin, fin] = greedyROI2d(Y_ds, nr, params);

% refine with HALS
fprintf('Refining initial estimates with HALS \n');
[Ain, Cin, bin, fin] = HALS_2d(Y_ds, full(Ain), Cin, bin, fin, params); 

%% upsample Ain, Cin, bin, fin
center = ssub*com(Ain,d1s,d2s); 
Ain = imresize(reshape(Ain, d1s, d2s, sum(nr)), [d1, d2]);
Ain = reshape(Ain, d1*d2, []); 
bin = imresize(reshape(bin, d1s, d2s), [d1, d2]);
bin = reshape(bin, [], 1); 
Cin = imresize(Cin, [sum(nr), Ts*tsub]);
fin = imresize(fin, [1, Ts*tsub]);
if T ~= Ts*tsub
    Cin = padarray(Cin, [0, T-Ts*tsub], 'post'); 
    fin = padarray(fin, [0, T-Ts*tsub], fin(end), 'post');   
end