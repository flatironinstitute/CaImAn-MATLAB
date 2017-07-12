clear;
addpath(genpath('../ca_source_extraction'));  % add packages to matlab path
addpath(genpath('../NoRMCorre'));
gcp;    % start a local cluster

filename = 'use_cases/CajalCourse2017/images_all.tif';
%% read file and determine dynamic range
Y = read_file(filename);
[d1,d2,T] = size(Y);    % dimensions of file
Y = Y - min(Y(:));      % remove negative offset

minY = quantile(Y(1:1e7),0.0005);
maxY = quantile(Y(1:1e7),1-0.0005);

%%  play raw movie 
play_movie({Y},{'Raw data'},minY,maxY)
%% perform motion correction (start with rigid)
tic
options_rg = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',100,'max_shift',15,'us_fac',50);
[M_rg,shifts_rg,template_rg] = normcorre_batch(Y,options_rg); 
time_mc_rig = toc
%% view data (press key to interrupt playing)
tsub = 5;   % downsampling factor
Y_sub = downsample_data(Y,'time',tsub);
M_rgs = downsample_data(M_rg,'time',tsub);

play_movie({Y_sub,M_rgs},{'Raw data (downsampled)','Rigidly corrected data (downsampled)'},minY,maxY)
%% clear memory


%% perform non-rigid motion correction    
%this does not improve the results
                
%% compute some metrics
[cY,mY,vY] = motion_metrics(Y,options_rg.max_shift);
[cM_rg,mM_rg,vM_rg] = motion_metrics(M_rg,options_rg.max_shift);
%% plot metrics
figure;
    ax(1) = subplot(2,2,1); imagesc(mY,[minY,maxY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax(2) = subplot(2,2,2); imagesc(mM_rg,[minY,maxY]); axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,2,3); plot(1:T,cY,1:T,cM_rg); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,2,4); scatter(cY,cM_rg); hold on; plot([0.9*min(cY),1.05*max(cM_rg)],[0.9*min(cY),1.05*max(cM_rg)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes(ax,'xy')
%%
clear Y Y_sub M_rgs
%% correlation image
Cn = correlation_image_max(M_rg);
imagesc(Cn)
axis image
%% now perform source extraction by splitting the FOV in patches

sizY = size(M_rg);
patch_size = [30,30];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 5;                                            % number of components to be found
tau = 5;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'cluster_pixels',false,...  
    'ssub',1,...                                % downsample factor in space
    'tsub',2,...                                % downsample factor in time
    'merge_thr',0.8,...                         % merging threshold
    'gSig',tau,... 
    'gnb',2,...
    'spatial_method','regularized'...
    );

%% run CNMF algorithm on patches and combine
tic;
[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(M_rg,K,patches,tau,p,options);
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(M_rg,A,C,b,f,YrA,options);
toc
%% restart pool here

%% a simple GUI
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
plot_components_GUI(M_nr,A(:,keep),C(keep,:),b,f,Cn,options);

%% refine temporal components
A_keep = A(:,keep);
C_keep = C(keep,:);
[C2,f2,P2,S2,YrA2] = update_temporal_components(reshape(M_rg,[],T),A_keep,b,C_keep,f,P,options);

git %% detrend fluorescence and extract DF/F values
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
figure;plot(1:T,F_dff(i,:),'-k'); hold all; plot(1:T,C_dec(i,:),'r','linewidth',2);
    spt = find(S(i,:));
    if spt(1) == 1; spt(1) = []; end
    hold on; scatter(spt,repmat(-0.25,1,length(spt)),'m*')
    title(['Component ',num2str(i)]);
    ls ..Cajal
    legend('Fluorescence DF/F','Deconvolved','Spikes')