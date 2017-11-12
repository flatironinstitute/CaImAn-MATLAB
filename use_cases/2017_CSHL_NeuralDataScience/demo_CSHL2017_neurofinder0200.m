clear all classes;
addpath(genpath('../../../ca_source_extraction'));  % add packages to matlab path
addpath(genpath('../../../NoRMCorre'));
gcp;    % start a local cluster

foldername = '/Users/epnevmatikakis/Documents/Ca_datasets/Neurofinder/neurofinder.02.00/images';
    % change foldername to where the data is saved 
files = subdir(fullfile(foldername,'image*.tif*'));
numFiles = length(files);

%% read neurofinder frames and combine them in a single tiff file
tic
if exist(fullfile(foldername,'neurofinder0200.tif'),'file')
    Ycon = read_file(fullfile(foldername,'neurofinder0200.tif'),1,1000);    
        % just read a few frames if it exists
else
    Ycon = read_neurofinder(foldername,fullfile(foldername,'neurofinder0200.tif'));
end
toc
%% get dynamic range for showing movies
minY = quantile(Ycon(1:1e7),0.002);
maxY = quantile(Ycon(1:1e7),1-0.002);
%%
play_movie({Ycon},{'Y'},minY,maxY)
%clear Ycon
%% perform rigid motion correction and save output as a tiff file
options_rg = NoRMCorreSetParms('d1',512,'d2',512,...
        'bin_width',200,'max_shift',15,'output_type','tif',...
        'tiff_filename',fullfile(foldername,'neurofinder0200_rig.tif'));

if ~exist(fullfile(foldername,'neurofinder0200_rig.tif'),'file') % check if file exists
    [M_rg,shifts_rg,template_rg] = normcorre_batch(fullfile(foldername,'neurofinder0200.tif'),options_rg); 
end


%% construct a memory mapped file

tic;
if ~exist(fullfile(foldername,'neurofinder0200_rig.mat'),'file')
    data = memmap_file(options_rg.tiff_filename);
else
    data = matfile(fullfile(foldername,'neurofinder0200_rig.mat'));
end
toc
%% now perform source extraction by splitting the FOV in patches

sizY = size(data,'Y');
patch_size = [42,42];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [5,5];                      % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 6;                                            % number of components to be found
tau = 5;                                          % std of gaussian kernel (half size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',2,...                                % downsample factor in space
    'tsub',4,...                                % downsample factor in time
    'merge_thr',0.8,...                         % merging threshold
    'gSig',tau,... 
    'gnb',2,...                                 % number of background components
    'spatial_method','regularized'...
    );


%% run CNMF algorithm on patches, combine results and classify components

tic;
[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);
toc;

%% compute max-correlation image (just for visualization purposes)
Cn = correlation_image_max(data);

%% a simple GUI for further classification
Coor = plot_contours(A,Cn,options,1,[]); close;

%% Evalute components to screen for false positive events
tic;
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(data,A,C,b,f,YrA,options);
toc

traces = prctfilt(C+YrA,8,1000,100);
ROIvars.fitness = compute_event_exceptionality(traces,0);
ROIvars.fitness_delta = compute_event_exceptionality(diff(traces,[],2),0);


%% view contour plots of selected and rejected components (optional)
keep = (fitness_delta < - 25 | fitness < - 25 | ROIvars.rval_space> 0.75) & ROIvars.sizeA > 10 & ROIvars.sizeA < 320  ; 
throw = ~keep;
figure;
    ax1 = subplot(121); plot_contours(A(:,keep),Cn,options,0,[],Coor,1,find(keep)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A(:,throw),Cn,options,0,[],Coor,1,find(throw));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')

%% run GUI for modifying component selection (optional, close twice to save values)
run_GUI = true;
if run_GUI
    %Coor = plot_contours(A,Cn,options,1); close;
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};    
end  

%% keep only the selected components and refine them
A_keep = A(:,keep);
C_keep = C(keep,:);
tic;
options.nb = size(b,2);
[A2,b2,C2,P2] = update_spatial_components(data,C_keep,f,[A_keep,b],P,options);
P2.p = 0;
[C2,f2,P2,S2,YrA2] = update_temporal_components(data,A2,b2,C2,f,P2,options);
toc

Coor2 = plot_contours(A2,Cn,options,1,[]); close
%% Evalute components to screen for false positive 

tic;
[ROIvars2.rval_space,ROIvars2.rval_time,ROIvars2.max_pr,ROIvars2.sizeA,keep2] = classify_components(data,A2,C2,b2,f2,YrA2,options);
toc
traces2 = prctfilt(C2+YrA2,8,1000,100);

fitness2 = compute_event_exceptionality(traces2,0);
fitness_delta2 = compute_event_exceptionality(diff(traces2,[],2),0);

%%
keep2 = (fitness_delta2 < - 30 | fitness2 < - 35 | ROIvars2.rval_space> 0.8) & ROIvars2.sizeA > 20 & ROIvars2.sizeA < 200 & ROIvars2.rval_time> 0.75 ; 
throw2 = ~keep2;
figure;
    ax1 = subplot(121); plot_contours(A2(:,keep2),Cn,options,0,[],Coor2,1,find(keep2)); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A2(:,throw2),Cn,options,0,[],Coor2,1,find(throw2));title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')

%% inspect components
A_keep2 = A2(:,keep2);
C_keep2 = C2(keep2,:);

plot_components_GUI(data,A_keep2,C_keep2,b2,f2,Cn,options);

%% detrend fluorescence and extract DF/F values

options.df_window = 1000; 
[F_dff,F0] = detrend_df_f(A_keep2,b2,C_keep2,f2,YrA2(keep2,:),options);

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