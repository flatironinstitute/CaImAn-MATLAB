clear;
%% load file

path_to_package = '';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
filename = '';      % path to stack tiff file
%foldername = '';    % path to folder that contains a sequence of tiff files

if exist([filename(1:end-3),'mat'],'file')
    data = matfile([filename(1:end-3),'mat'],'Writable',true);
else
    sframe=1;						% user input: first frame to read (optional, default 1)
    num2read=[];					% user input: how many frames to read   (optional, default until the end)
    chunksize=1000;                 % user input: read and map input in chunks (optional, default read all at once)
    data = memmap_file(filename,sframe,num2read,chunksize);
    %data = memmap_file_sequence(foldername);
end

%% Set parameters
sizY = size(data,'Y');                  % size of data matrix
patch_size = [32,32];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [4,4];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 8;                                            % number of components to be found
tau = 8;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
sizY = data.sizY;

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',2,...
    'tsub',4,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% order and plot

[A_or,C_or,S_or,P] = order_ROIs(A,C,S,P); % order components

contour_threshold = 0.95;                 % amount of energy used for each component to construct contour plot
figure;
[Coor,json_file] = plot_contours(A_or,reshape(P.sn,sizY(1),sizY(2)),contour_threshold,1); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)