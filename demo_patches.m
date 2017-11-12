% demo script for splitting the field of view in patches and processing in parallel
% with or without memory mapping. See also run_pipeline.m for the complete  
% pre-processing pipeline of large datasets

clear;
%% setup path to file and package

path_to_package = '../ca_source_extraction';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
filename = '/Users/epnevmatikakis/Documents/Ca_datasets/Neurofinder/neurofinder.02.00/images/neurofinder0200_rig.tif';      
        % path to file (assumed motion corrected)
        
is_memmaped = false;        % choose whether you want to load the file in memory or not

%% load file

if is_memmaped
    if exist([filename(1:end-3),'mat'],'file')
        data = matfile([filename(1:end-3),'mat'],'Writable',true);
    else
        sframe=1;						% user input: first frame to read (optional, default 1)
        num2read=[];					% user input: how many frames to read   (optional, default until the end)
        chunksize=5000;                 % user input: read and map input in chunks (optional, default read all at once)
        data = memmap_file(filename,sframe,num2read,chunksize);
        %data = memmap_file_sequence(foldername);
    end
    sizY = size(data,'Y');                    % size of data matrix
else
    T = 2000;                                 % load only a part of the file due to memory reasons
    data = read_file(filename,1,T);
    sizY = size(data);
end
    
%% Set parameters
patch_size = [32,32];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [6,6];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 10;                  % number of components to be found
tau = 7;                 % std of gaussian kernel (size of neuron) 
p = 2;                   % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;         % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','dilate',...                % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'ssub',2,...
    'tsub',1,...
    'p',p,...                                   % order of AR dynamics
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'spatial_method','regularized',...
    'min_SNR',2);

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% classify components 

[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA] = classify_components(data,A,C,b,f,YrA,options);
ind_corr = ROIvars.rval_space > options.space_thresh;  % space correlation

try  % matlab 2017b or later is needed for the CNN classifier
    [ind_cnn,value] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);
end

N_samples_exc = ceil(options.fr*options.decay_time);  % event exceptionality
fitness = compute_event_exceptionality(C+YrA,N_samples_exc,options.robust_std);
ind_exc = (fitness < log(normcdf(-options.min_SNR))*N_samples_exc);

keep = (ind_corr & ind_cnn) | ind_exc;

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
A_throw = A(:,~keep);
C_throw = C(~keep,:);
A_keep = A(:,keep);
C_keep = C(keep,:);
options.p = 2;      % perform deconvolution
P.p = 2;
[C2,f2,P2,S2,YrA2] = update_temporal_components_fast(data,A_keep,b,C_keep,f,P,options);

%% plot results
options.sx = 64;
plot_components_GUI(double(data),A,C,b,f,Cn,options);
