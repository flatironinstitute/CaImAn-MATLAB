clear;
% same demo as demo_script.m but using the clas @Sources2D
%% load file

addpath(genpath('utilities'));
addpath(genpath('../cvx'));
addpath(genpath('../constrained_foopsi'));


%nam = 'demoMovie.tif';          % insert path to tiff stack here
data_root='Q:\data\2photon\reg';
data_folder='141221_DM044_2P_DM\run02_blank\';
nam = fullfile(data_root,data_folder,'run02_blank_green_reg000001.tif');

switch getenv('computername')
    case 'BKRUNCH'
        
        options.req_NumWorkers=8;
        options.nBlocks=64;
        options.frame_rate_actual=30;
        options.frame_rate_req=4; % Hz: defines resampling rate => 4Hz = look at correlated activity in 250 ms bins
        red=[linspace(0,1,256)' zeros(256,1) zeros(256,1)];
        green=[zeros(256,1) linspace(0,1,256)' zeros(256,1)];

        if isempty(gcp('nocreate'))
            PP=parpool(options.req_NumWorkers);
        end
        
        exp_name=data_folder;
        expt = frGetExpt(exp_name);
        
        filelist = dir(fullfile(expt.dirs.reggreenpn,'*.tif'));
        nFiles = length(filelist);
        
        options.nBlocks=min([options.nBlocks nFiles]);
        stack = readtiff(expt.dirs.regrootpn,1:options.nBlocks);
        
        %% Reduce frame_rate
        %nFrames_req=32*options.frame_rate_req; % Hz
        %bin_frames=round(size(stack,3)/nFrames_req)
        bin_frames=round(options.frame_rate_actual/options.frame_rate_req)
        dec = stackGroupProject(stack,bin_frames,'sum');
        Y = double(dec);
    otherwise
        sframe=1;						% user input: first frame to read (optional, default 1)
        num2read=2000;					% user input: how many frames to read   (optional, default until the end)
        Y = bigread2(nam,sframe,num2read);
end

if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% Set parameters

K = 200;                                          % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron)
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold

obj = Sources2D;
updateParams(obj,...
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr...                    % merging threshold
    );

%% Data pre-processing
Y = preprocess(obj,Y,p);

%% fast initialization of spatial components using greedyROI and HALS
center = initComponents(obj, Y, K, tau);

% display centers of found components
Cn =  correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
axis equal; axis tight; hold all;
scatter(center(:,2),center(:,1),'mo');
title('Center of ROIs found from initialization algorithm');
drawnow;

%% update spatial components
Yr = reshape(Y,d,T);
clear Y;
updateSpatial(obj, Yr);

%% update temporal components
Y_res = updateTemporal(obj, Yr);

%% merge found components
Apr = obj.A;    % store non-merged components
Cpr = obj.C;
[K_m, merged_ROIs] = merge(obj, Y_res);
display_merging = 1; % flag for displaying merging example
if display_merging
    i = 1; randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
    set(gcf,'Position',[300,300,(ln+2)*300,300]);
    for j = 1:ln
        subplot(1,ln+2,j); imagesc(reshape(Apr(:,merged_ROIs{i}(j)),d1,d2));
        title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
    end
    subplot(1,ln+2,ln+1); imagesc(reshape(obj.A(:,K_m-length(merged_ROIs)+i),d1,d2));
    title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight;
    subplot(1,ln+2,ln+2);
    plot(1:T,(diag(max(Cpr(merged_ROIs{i},:),[],2))\Cpr(merged_ROIs{i},:))');
    hold all; plot(1:T,obj.C(K_m-length(merged_ROIs)+i,:)/max(obj.C(K_m-length(merged_ROIs)+i,:)),'--k')
    title('Temporal Components','fontsize',16,'fontweight','bold')
    drawnow;
end

%% repeat
updateSpatial(obj, Yr);
Y_res = updateTemporal(obj, Yr);
[C_df, ~, S_df] = extractDF_F(obj, Yr, K_m+1);

%% do some plotting
[srt] = orderROIs(obj);     % order components
contour_threshold = 0.95;   % amount of energy used for each comp vonent to construct contour plot
if 0
    figure;
    [json_file] = viewContours(obj, Cn, contour_threshold, 1);
    pause;
    %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
    plotComponents(obj, Yr, Cn);     % display all components
    %% make movie
    makePatchVideo(obj, Yr)
end

%%
figure(55)
imagesc(C_df)
colormap(green)