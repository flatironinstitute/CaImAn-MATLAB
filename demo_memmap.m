% demo script for splitting the field of view in patches and processing in parallel
% through memory mapping. See also run_pipeline.m for the complete  
% pre-processing pipeline of large datasets

clear;
%% load file

path_to_package = '../ca_source_extraction';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
filename = '/Users/epnevmatikakis/Documents/Ca_datasets/N5TailPuff_0_0___0,000_0,000_1_1_m_r_nb.tif';      % path to stack tiff file
%foldername = '';    % path to folder that contains a sequence of tiff files
%%
if exist([filename(1:end-3),'mat'],'file')
    data = matfile([filename(1:end-3),'mat'],'Writable',true);
else
    sframe=1;						% user input: first frame to read (optional, default 1)
    num2read=[];					% user input: how many frames to read   (optional, default until the end)
    chunksize=5000;                 % user input: read and map input in chunks (optional, default read all at once)
    data = memmap_file(filename,sframe,num2read,chunksize);
    %data = memmap_file_sequence(foldername);
end

%% Set parameters
sizY = size(Y);                  % size of data matrix
patch_size = [445,582];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [32,32];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 100;                                            % number of components to be found
tau = 16;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','dilate','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'init_method','HALS',...
    'cluster_pixels',false,...
    'ssub',2,...
    'tsub',1,...
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'spatial_method','regularized'...
    );

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(Y,K,patches,tau,p,options);

%% classify components
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);

%% run GUI for modifying component selection (optional, close twice to save values)
Cn = reshape(P.sn,sizY(1),sizY(2));  % background image for plotting
run_GUI = false;
if run_GUI
    Coor = plot_contours(A,Cn,options,1); close;
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};    
end

%% re-estimate temporal components

A_keep = A(:,keep);
C_keep = C(keep,:);
options.p = 2;      % perform deconvolution
[C2,f2,P2,S2,YrA2] = update_temporal_components_fast(data,A_keep,b,C_keep,f,P,options);

%% plot results
options.sx = 64;
plot_components_GUI(double(Y),A,C,b,f,Cn,options);