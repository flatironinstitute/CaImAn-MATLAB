clear;
addpath(genpath('../../../ca_source_extraction'));  % add packages to matlab path
addpath(genpath('../../../NoRMCorre'));
gcp;    % start a local cluster

foldername = '/Users/epnevmatikakis/Documents/Ca_datasets/Neurofinder/neurofinder.02.00/images';
    % change foldername to where the data is saved 
files = subdir(fullfile(foldername,'*.tif*'));
numFiles = length(files);

%% concatenate files (this will take some time)
tic;
Ycon = concatenate_files(files,fullfile(foldername,'neurofinder0200.tif'),'tif');
toc
%% construct a memory mapped file
tic;
data = memmap_file(fullfile(foldername,'neurofinder0200.tif'));
toc

%% now perform source extraction by splitting the FOV in patches

sizY = size(data,'Y');
patch_size = [48,48];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 6;                                            % number of components to be found
tau = 7;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'cluster_pixels',false,...  
    'ssub',2,...                                % downsample factor in space
    'tsub',4,...                                % downsample factor in time
    'merge_thr',0.8,...                         % merging threshold
    'gSig',tau,... 
    'gnb',2,...
    'spatial_method','regularized'...
    );

%% run CNMF algorithm on patches, combine results and classify components
tic;
[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(data,A,C,b,f,YrA,options);
toc

%% a simple GUI for further classification
Coor = plot_contours(A,Cn,options,1); close;
% run_GUI = 0;
% if run_GUI
%     GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
%     options = GUIout{2};
%     keep = GUIout{3};    
% end

%% view contour plots of selected and rejected components (optional)
keep = (ROIvars.rval_space>.7 & ROIvars.rval_time>0);
throw = ~keep;
figure;
    ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')
%% inspect components
plot_components_GUI(data,A(:,keep),C(keep,:),b,f,Cn,options);

%% refine temporal components
A_keep = A(:,keep);
C_keep = C(keep,:);
[C2,f2,P2,S2,YrA2] = update_temporal_components(data,A_keep,b,C_keep,f,P,options);

%% detrend fluorescence and extract DF/F values
df_percentile = 30;
window = 1000; 

F = diag(sum(A_keep.^2))*(C2 + YrA2);  % fluorescence
Fd = prctfilt(F,df_percentile,window);                      % detrended fluorescence
Bc = prctfilt((A_keep'*b)*f2,30,1000,300,0) + (F-Fd);       % background + baseline for each component
F_dff = Fd./Bc;

%% deconvolve data

nNeurons = size(F_dff,1);
C_dec = zeros(size(F_dff));
S = zeros(size(F_dff));
kernels = cell(nNeurons,1);
min_sp = 3;    % find spikes resulting in transients above min_sp x noise level

for i = 1:nNeurons
    [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
end

%% plot a random component
i = randi(nNeurons);
T = sizY(end);
figure;plot(1:T,F_dff(i,:),'-k'); hold all; plot(1:T,C_dec(i,:),'r','linewidth',2);
    spt = find(S(i,:));
    if spt(1) == 1; spt(1) = []; end
    hold on; scatter(spt,repmat(-0.25,1,length(spt)),'m*')
    title(['Component ',num2str(i)]);
    legend('Fluorescence DF/F','Deconvolved','Spikes')