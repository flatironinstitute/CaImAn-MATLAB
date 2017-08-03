clear;
gcp;    % start a local cluster

filename = 'demoSue2x.tif';
if ~exist(filename,'file');
    url = 'https://www.dropbox.com/s/36xdfd28eone0hj/demoSue2x.tif?dl=1';
    fprintf('downloading the file...');
    outfilename = websave(filename,url);
    fprintf('done. \n');
end

addpath(genpath('../../../ca_source_extraction-master'));  % add packages to matlab path
addpath(genpath('../../../NoRMCorre-master'));
addpath(genpath('../../../ca_source_extraction'));  % add packages to matlab path
addpath(genpath('../../../NoRMCorre'));

%% read file and determine dynamic range

Y = read_file(filename);
[d1,d2,T] = size(Y);    % dimensions of file
Y = Y - min(Y(:));      % remove negative offset

minY = quantile(Y(1:1e7),0.0005);
maxY = quantile(Y(1:1e7),1-0.0005);

%%  view data
figure;play_movie({Y},{'raw data'},minY,maxY);

    
%% perform motion correction (start with rigid)
% parameters motion correction
% 'd1','d2': size of FOV
% 'bin_width': how often to update the template
% 'max_shift': maximum allowed rigid shift

options_rg = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',100,'max_shift',15);

[M_rg,shifts_rg,template_rg] = normcorre_batch(Y,options_rg); 

%% view data
tsub = 5;   % downsampling factor (only for display purposes)
Y_sub = downsample_data(Y,'time',tsub);
M_rgs = downsample_data(M_rg,'time',tsub);

play_movie({Y_sub,M_rgs},{'raw data','rigid'},minY,maxY);

%% perform non-rigid motion correction    
% parameters motion correction
% 'd1','d2': size FOV movie
% 'grid_size','overlap_pre': parameters regulating size of patch (size patch ~ (grid_size + 2*overlap_pre))
% 'mot_uf': upsampling factor of the grid for shift application
% 'bin_width': how often to update the template
% 'max_shift': maximum allowed rigid shift
% 'max_dev': maximum deviation allowed for each patch from the rigid shift value

options_nr = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),...
                    'grid_size',[48,48],'mot_uf',4,'overlap_pre',[16,16],...
                    'bin_width',100,'max_shift',15,'max_dev',8);

[M_nr,shifts_nr,template_nr] = normcorre_batch(Y,options_nr,template_rg); 

%% view (downsampled) data

M_nrs = downsample_data(M_nr,'time',tsub);
play_movie({Y_sub,M_rgs,M_nrs},{'raw data','rigid','pw-rigid'},minY,maxY);

%% compute some metrics for motion correction quality assessment

[cY,mY,vY] = motion_metrics(Y,options_rg.max_shift);
[cM_rg,mM_rg,vM_rg] = motion_metrics(M_rg,options_rg.max_shift);
[cM_nr,mM_nr,vM_nr] = motion_metrics(M_nr,options_rg.max_shift);

%% plot shifts        

shifts_r = squeeze(cat(3,shifts_rg(:).shifts));
shifts_n = cat(ndims(shifts_nr(1).shifts)+1,shifts_nr(:).shifts);
shifts_n = reshape(shifts_n,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_n(:,2,:))';
shifts_y = squeeze(shifts_n(:,1,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM_rg,1:T,cM_nr); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')                


%% plot metrics
figure;
    ax(1) = subplot(2,3,1); imagesc(mY,[minY,maxY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax(2) = subplot(2,3,2); imagesc(mM_rg,[minY,maxY]); axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    ax(3) = subplot(2,3,3); imagesc(mM_nr,[minY,maxY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM_rg,1:T,cM_nr); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM_rg); hold on; plot([0.9*min(cY),1.05*max(cM_rg)],[0.9*min(cY),1.05*max(cM_rg)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM_rg,cM_nr); hold on; plot([0.95*min(cM_rg),1.05*max(cM_nr)],[0.95*min(cM_rg),1.05*max(cM_nr)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes(ax,'xy')


%% now perform source extraction by splitting the FOV in patches

sizY = size(M_nr);
patch_size = [30,30];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [8,8];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 4;                                            % number of components to be found
tau = 4;                                          % std of gaussian kernel (half size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'ssub',1,...                                % downsample in space
    'tsub',2,...                                % downsample in time
    'merge_thr',0.8,...                         % merging threshold
    'gSig',tau,... 
    'gnb',2,...                                 % number of background components
    'spatial_method','regularized'...
    );

%% run CNMF algorithm on patches and combine
tic;
[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(M_nr,K,patches,tau,p,options);
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(M_nr,A,C,b,f,YrA,options);
toc

%% a simple GUI
Cn = correlation_image_max(M_nr);
Coor = plot_contours(A,Cn,options,1); close;
run_GUI = false;
if run_GUI
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};    
end

%% view contour plots of selected and rejected components (optional)
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
P.p = 0;
[C2,f2,P2,S2,YrA2] = update_temporal_components(reshape(M_nr,[],T),A_keep,b,C_keep,f,P,options);

%% detrend fluorescence and extract DF/F values

options.df_window = 1000; 
[F_dff,F0] = detrend_df_f(A_keep,b,C2,f2,YrA2,options);

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

figure;plot(1:T,F_dff(i,:),'--k'); hold all; plot(1:T,C_dec(i,:),'r','linewidth',2);
    spt = find(S(i,:));
    if spt(1) == 1; spt(1) = []; end
    hold on; scatter(spt,repmat(-0.25,1,length(spt)),'m*')
    title(['Component ',num2str(i)]);
    
    legend('Fluorescence DF/F','Deconvolved','Spikes')