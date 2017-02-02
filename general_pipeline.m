% TAXIDIS pipeline
clear;
addpath(genpath('../NoRMCorre'));               % add the NoRMCorre motion correction package to MATLAB path
gcp;        % start a parallel engine
foldername = '/mnt/ceph/users/epnevmatikakis/Ca_datasets/Taxidis/mouse1-multi-sessions';   
        % folder where all the files are located. Currently supported .tif,
        % .hdf5, .raw, .avi, and .mat files
files = subdir(fullfile(foldername,'*.raw'));   % list of filenames (will search all subdirectories)
FOV = [512,512];
numFiles = length(files);

%% motion correct (and save registered h5 files as 2d matrices (to be used in the end)..)
% register files one by one. use template obtained from file n to
% initialize template of file n + 1; 

non_rigid = true; % flag for non-rigid motion correction

template = [];
for i = 1:numFiles
    name = files(i).name;
    if non_rigid
        options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512,'grid_size',[128,128],...
            'overlap_pre',64,'mot_uf',4,'bin_width',100,'max_shift',24,'max_dev',8,'us_fac',50,...
            'output_type','h5','h5_filename',[name(1:end-4),'_nr.h5']);
        [M,shifts,template] = normcorre_batch(name,options_rigid,template); 
        save([name(1:end-4),'_shifts_nr.mat'],'shifts','-v7.3');           % save shifts of each file at the respective subfolder
    else    % perform rigid motion correction (faster, could be less accurate)
        options_rigid = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',100,'max_shift',32,...
            'output_type','h5','h5_filename',[name(1:end-4),'_rig.h5']);
        [M,shifts,template] = normcorre_batch(name,options_rigid,template); 
        save([name(1:end-4),'_shifts_rig.mat'],'shifts','-v7.3');           % save shifts of each file at the respective subfolder
    end
end

%% downsample h5 files and save into a single memory mapped matlab file

if non_rigid
    h5_files = subdir(fullfile(foldername,'*_nr.h5'));  % list of h5 files (modify 
else
    h5_files = subdir(fullfile(foldername,'*_rig.h5'));
end

tsub = 5;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
ds_filename = [foldername,'/ds_data.mat'];
data = matfile(ds_filename,'Writable',true);
data.Y  = zeros([FOV,0],'uint16');
data.Yr = zeros([prod(FOV),0],'uint16');
data.sizY = [FOV,0];

batch_size = 2000;                               % read chunks of that size
batch_size = round(batch_size/tsub)*tsub;        % make sure batch_size is divisble by tsub
Ts = zeros(numFiles,1);                          % store length of each file
cnt = 0;                                         % number of frames processed so far
tt1 = tic;
for i = 1:numFiles
    name = h5_files(i).name;
    info = h5info(name);
    dims = info.Datasets.Dataspace.Size;
    ndimsY = length(dims);                       % number of dimensions (data array might be already reshaped)
    Ts(i) = dims(end);
    Ysub = zeros(FOV(1),FOV(2),floor(Ts(i)/tsub),'uint16');
    data.Y(FOV(1),FOV(2),sum(floor(Ts/tsub))) = uint16(0);
    data.Yr(prod(FOV),sum(floor(Ts/tsub))) = uint16(0);
    cnt_sub = 0;
    for t = 1:batch_size:Ts(i)
        Y = bigread2(name,t,min(batch_size,Ts(i)-t+1));                
        ln = size(Y,ndimsY);
        Y = reshape(Y,[FOV,ln]);
        Y = uint16(downsample_data(uint16(Y),'time',tsub));
        ln = size(Y,3);
        Ysub(:,:,cnt_sub+1:cnt_sub+ln) = Y;
        cnt_sub = cnt_sub + ln;
    end
    data.Y(:,:,cnt+1:cnt+cnt_sub) = Ysub;
    data.Yr(:,cnt+1:cnt+cnt_sub) = reshape(Ysub,[],cnt_sub);
    toc(tt1);
    cnt = cnt + cnt_sub;
    data.sizY(1,3) = cnt;
end

%% now run CNMF on patches on the downsampled file, set parameters first

sizY = data.sizY;                       % size of data matrix
patch_size = [40,40];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 6;                                            % number of components to be found
tau = 8;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
sizY = data.sizY;

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'cluster_pixels',false,...
    'ssub',2,...                                % spatial downsampling when processing
    'tsub',4,...                                % further temporal downsampling when processing
    'fudge_factor',0.96,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'max_size_thr',300,'min_size_thr',10,...    % max/min acceptable size for each component
    'spatial_method','constrained',...
    'df_prctile',50,...                         % take the median of background fluorescence to compute baseline fluorescence 
    'fr',30/tsub...
    );

%% Run on patches (around 15 minutes)

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% compute correlation image on a small sample of the data (optional - for visualization purposes) 
Cn = correlation_image(single(data.Y(:,:,1:2000)),8);
%% get contour plots for all components
CC = plot_contours(A,Cn,options,0); close;
%% classify components (around 3 minutes)
options.space_thresh = 0.5;
options.time_thresh = 0.5;
[rval_space,rval_time,max_pr,sizeA,keep] = classify_components(data,A,C,b,f,YrA,options);

%% modify thresholds to examine which parameters are kept

keep = (rval_space > 0.3) & (rval_time > options.time_thresh) & (max_pr > 0.9) & (sizeA >= options.min_size_thr) & (sizeA <= options.max_size_thr);
    % play with thresholds to modify the components being selected
%% view contour plots of selected and rejected components (optional)

throw = ~keep;
figure;
    ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],CC,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],CC,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')
%%

gui_components(rval_space,rval_time,max_pr,sizeA,A,Cn,options)
%% keep only the active components    
A_keep = A(:,keep);
C_keep = C(keep,:);

%% deconvolve (downsampled) temporal components plot GUI with components (optional)

%tic;
%[C_keep,f_keep,Pk,Sk,YrAk] = update_temporal_components_fast(data,A_keep,b,C_keep,f,P,options);
%toc

plot_components_GUI(data,A_keep,C_keep,b,f,Cn,options)

%% extract fluorescence and DF/F on native temporal resolution (13-14 minutes)
% C is deconvolved activity, C + YrA is non-deconvolved fluorescence 
% F_df is the DF/F computed on the non-deconvolved fluorescence

P.p = 0;                    % order of dynamics. Set P.p = 0 for no deconvolution at the moment
C_us = cell(numFiles,1);    % cell array for thresholded fluorescence
f_us = cell(numFiles,1);    % cell array for temporal background
P_us = cell(numFiles,1);  
S_us = cell(numFiles,1);
YrA_us = cell(numFiles,1);  % 
b_us = cell(numFiles,1);    % cell array for spatial background
F_df = cell(numFiles,1);    % cell array for DF/F values
Df = cell(numFiles,1);
tt1 = tic;   
for i = 1:numFiles    
    int = sum(floor(Ts(1:i-1)/tsub))+1:sum(floor(Ts(1:i)/tsub));
    Cin = imresize([C_keep(:,int);f(:,int)],[size(C_keep,1)+size(f,1),Ts(i)]);
    [C_us{i},f_us{i},P_us{i},S_us{i},YrA_us{i}] = update_temporal_components_fast(h5_files(i).name,A_keep,b,Cin(1:end-1,:),Cin(end,:),P,options);
    b_us{i} = max(mm_fun(f_us{i},h5_files(i).name) - A_keep*(C_us{i}*f_us{i}'),0)/norm(f_us{i})^2;
    [F_df{i},Df{i}] = extract_DF_F_new(A_keep,C_us{i}+YrA_us{i},b_us{i},f_us{i},P_us{i},options);
    toc(tt1);
end
F_us = cellfun(@plus,C_us,YrA_us,'un',0);       % cell array for 

%% perform deconvolution