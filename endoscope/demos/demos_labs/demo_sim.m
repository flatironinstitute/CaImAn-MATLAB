%% clear workspace
clear; clc; close all;

%% select data and map it with the memory
if ~exist('nam', 'var') || isempty(nam)
    try
        load .dir.mat; %load previous path
    catch
        dir_nm = [cd(), filesep]; %use the current path
    end
    [file_nm, dir_nm] = uigetfile(fullfile(dir_nm, '*.tif;*.mat'));
    if dir_nm~=0
        save .dir.mat dir_nm;
    else
        fprintf('no file was selecteD. STOP!\N');
        return;
    end
    nam = [dir_nm, file_nm];  % full name of the data file
    [~, file_nm, file_type] = fileparts(nam);
end

% convert the data to mat file
nam_mat = [dir_nm, file_nm, '.mat']; 
if strcmpi(file_type, '.mat')
    fprintf('The selected file is *.mat file\n'); 
elseif  exist(nam_mat', 'file')
    % the selected file has been converted to *.mat file already 
    fprintf('The selected file has been replaced with its *.mat version\n'); 
elseif or(strcmpi(file_type, '.tif'), strcmpi(file_type, '.tiff'))
    % convert 
    tic;
    fprintf('converting the selected file to *.mat version...\n'); 
    nam_mat = tif2mat(nam);
    fprintf('Time cost in converting data to *.mat file:     %.2f seconds\n', toc);
else
    fprintf('The selected file type was not supported yet! email me to get support (zhoupc1988@gmail.com)\n');
    return; 
end

data = matfile(nam_mat); 
Ysiz = data.Ysiz;
d1 = Ysiz(1);   %height
d2 = Ysiz(2);   %width
numFrame = Ysiz(3);    %total number of frames

fprintf('\nThe data has been mapped to RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1, d2, numFrame, prod(Ysiz)*8/(2^30));

%% create Source2D class object for storing results and parameters
neuron_raw = Sources2D('d1',d1,'d2',d2);   % dimensions of datasets
neuron_raw.Fs = 20;         % frame rate
neuron_raw.options.Fs = 20; 
ssub = 1;           % spatial downsampling factor 
tsub = 1;           % temporal downsampling factor 
neuron_raw.updateParams('ssub', ssub,...  % spatial downsampling factor
    'tsub', tsub, ...  %temporal downsampling factor
    'gSig', 4,... %width of the gaussian kernel, which can approximates the average neuron shape
    'gSiz', 15, ...% average size of a neuron
    'dist', 2, ... % maximum size of the neuron: dist*gSiz
    'search_method', 'ellipse', ... % searching method
    'merge_thr', 0.7, ... % threshold for merging neurons
    'bas_nonneg', 1);   % 1: positive baseline of each calcium traces; 0: any baseline

%% downsample data for fast and better initialization 
sframe=1;						% user input: first frame to read (optional, default:1)
num2read= numFrame;             % user input: how many frames to read   (optional, default: until the end)

tic; 
if and(ssub==1, tsub==1)
    neuron = neuron_raw; 
    Y = double(data.Y);  
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
else 
    [Y, neuron] = neuron_raw.load_data(nam_mat, sframe, num2read);
    [d1s,d2s, T] = size(Y);
    fprintf('\nThe data has been downsampled and loaded into RAM. It has %d X %d pixels X %d frames. \nLoading all data requires %.2f GB RAM\n\n', d1s, d2s, T, d1s*d2s*T*8/(2^30));
end 
Y = neuron.reshape(Y, 1);
neuron_raw.P.p = 2;      %order of AR model

fprintf('Time cost in downsapling data:     %.2f seconds\n', toc);

%% initialization of A, C
tic;
debug_on = true;        % debug mode 
save_avi = true;
neuron.options.min_corr = 0.85;
neuron.options.nk = 1; %round(T/(60*neuron.Fs)); % number of knots for spline basis, the interval between knots is 180 seconds
patch_par = 1; %[3, 3];  % divide the optical field into 3 X 3 patches and do initialization patch by patch
K = 100; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
neuron.options.bd = []; % boundaries to be removed due to motion correction 
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 
figure; 
imagesc(Cn); 
hold on; plot(center(:, 2), center(:, 1), 'or'); 
colormap; axis off tight equal; 

[~, srt] = sort(max(neuron.C, [], 2)./get_noise_fft(neuron.C), 'descend'); 
neuron.orderROIs(srt); 
neuron_init = neuron.copy(); 

%% merge neurons, order neurons and delete some low quality neurons (cell 0, before running iterative udpates)
[Ain, Cin] = neuron.snapshot();   % keep the initialization results
neuron.options.merge_thr = 0.5;     % threshold for merging neurons
% merged_ROI = neuron.quickMerge('C');  %merge neurons based on their temporal correlation (calcium traces)
merged_ROI = neuron.quickMerge('S');  %merge neurons based on their
% temporal correlation (spike count)
display_merge = true; 
if display_merge && ~isempty(merged_ROI)
    figure; 
    for m=1:length(merged_ROI)
        subplot(221);
        neuron.image(sum(Ain(:, merged_ROI{m}), 2));
        axis equal off tight; 
        subplot(2,2,3:4);
        plot(Cin(merged_ROI{m}, :)');
        axis tight;
        pause;
    end
end

[Cpnr, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend'); 
neuron.orderROIs(srt); 
[Ain, Cin] = neuron.snapshot();   % keep the initialization results
neuron.viewNeurons(); 

%% display contours of the neurons 
figure; 
neuron.viewContours(Cn, 0.8, 0); 
colormap winter; 
axis equal; axis off;
title('contours of estimated neurons');

%% udpate background (cell 1, the following three cells can be run iteratively)
% start approximating theb background
tic;
Ybg = Y-neuron.A*neuron.C; 
ssub = 2;   % downsample the data to improve the speed 
rr = round(neuron.options.gSiz*1.5);  % average neuron size, it will determine the neighbors for regressing each pixel's trace
active_px = [];% (sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
Ybg = neuron.localBG(Ybg, ssub, rr, active_px); % estiamte local background.
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc); 

% subtract the background from the raw data. 
Ysignal = Y - Ybg;
if ~isfield(neuron.P, 'sn')
    neuron.preprocess(Ysignal, 2);
end
% neuron.playVideo(Ysignal); % play the video data after subtracting the background components. 
%% update spatial components (cell 2), we can iteratively run cell 2& 3 for few times and then run cell 1
spatial_method = 'hals'; % methods for updating spatial components {'hals', 'lars'}, 
% hals is fast, there is no sparse constraint; lars is very slow, but the
% results are sparse. LARS is very slow. so run LARS in the last step 
tic;
% update spatial components with model Y = A*C
if strcmpi(spatial_method, 'lars')
    if ~isfield(neuron.P, 'sn')
        neuron.preprocess(Ysignal, 2);
    end
    neuron.updateSpatial_nb(Ysignal);
else
    IND = determine_search_location(neuron.A, 'dilate', neuron.options);
    neuron.A = HALS_spatial(Ysignal, neuron.A, neuron.C, IND, 10);
    %     neuron.post_process_spatial(); % uncomment this line to postprocess
    ind = find(sum(neuron.A, 1)==0);
    neuron.delete(ind);
    %     the results
end
fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);

%% update C  (cell 3)
temporal_method = 'hals'; % methods for updating temporal components {'hals', 'deconv'}
tic;
C0 = neuron.C;
if strcmpi(temporal_method, 'deconv')
    neuron.C = HALS_temporal(Ysignal, neuron.A, neuron.C, 10); 
    neuron.updateTemporal_nb(Ysignal); 
else
    neuron.C = HALS_temporal(Ysignal, neuron.A, neuron.C, 10);
end

fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);

%% pick neurons from the residual (cell 3.5). It's not alway necessary 
neuron.manual_add(Ysignal-neuron.A*neuron.C); 

%% Event detection (This should be the last step when you are satisfied with the neuron segmentation and the background estimation)
% algorithm: find local maximum within a window of 1 seconds first; if the
% difference between the local max and nearby minimum value is larger than
% thr_factor*noise, then label this local max as an event
thr_factor = 10; % the minimum calcium transients should be larger than thr_factor*noise
E = neuron.event_detection(thr_factor, neuron.Fs); 
neuron.viewNeurons([], 1.5*E); 

%% save results 
result_nm = [dir_nm, file_nm, '_results.mat']; 
neuron.save_results(result_nm); %save variable 'neuron' only. 
% neuron.save_results(result_nm, neuron.reshape(Ybg, 2)); % save background as well 

%% display neurons
figure;
neuron.viewNeurons();

%% display contours of the neurons 
figure; 
% neuron.viewContours(Cn, 0.9, 0);  % correlation image computed with
% spatially filtered data 
Cnn = correlation_image(Ysignal(:, 1:5:end), 4, d1, d2); 
neuron.viewContours(Cnn, 0.9, 0); % correlation image computed with background-subtracted data
colormap winter; 
axis equal; axis off;
title('contours of estimated neurons');
% plot contours with IDs
figure; 
plot_contours(neuron.A, Cn, 0.9, 1, [], neuron.Coor); 

%% check results by visualizing all spatial and temporal components one by one
folder_nm = [];%'neurons';
neuron.viewNeurons([], C, folder_nm);

%% check spatial and temporal components by playing movies
save_avi = false;
avi_name = 'play_movie.avi';
neuron.Cn = Cnn; 
neuron.runMovie(Ysignal, [0, 50], save_avi, avi_name);

%%

Yac = neuron.reshape(neuron.A*neuron.C, 2);
Ybg = neuron.reshape(Ybg, 2);
% Y = neuron.reshape(Y, 2);
Y = neuron.reshape(Y, 2);
ctr = round( neuron.estCenter());
% figure;
% neuron.viewContours(Cn, .5, 0);
% cell_IDs = [];

figure('position', [0,0, 1248, 600]);
% avi_file = VideoWriter('~/Dropbox/Public/Garret/day1/residual.avi');
save_avi = false;
if save_avi
    avi_file = VideoWriter([dir_nm, file_nm, '_results1.avi']);
    avi_file.open();
end
Ymax = 2000; %quantile(Y(1:1000:(d1*d2*T)), 0.9999);
ACmax = 50; %quantile(Yac(1:1000:(d1*d2*T)), 0.9999);

%     subplot(4,6, [5,6,11,12]); 
for m=1:5:T
    subplot(4,6, [1,2, 7, 8]);
    imagesc(Y(:, :,m), [0, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw data'); hold on; colorbar;
    
    subplot(4, 6, [1,2, 7, 8]+12);
         imagesc(Ybg(:, :, m), [0, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('background');
    colorbar;
    
    subplot(4,6, [3,4, 9,10]);
    imagesc(Y(:, :, m)-Ybg(:, :, m), [0, ACmax]); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw-background'); hold on; colorbar;
    
            
    subplot(4, 6, [3,4, 9,10]+12); 
    imagesc(Y(:, :, m)-Ybg(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('residual'); hold on; colorbar;

    subplot(4,6, [5,6,11,12]); 
        imagesc(neuron.reshape(neuron.A*neuron.C(:, m), 2), [0, ACmax]);
        %     imagesc(Ybg(:, :, m), [-50, 50]);
        text(d2/2-30, 10, sprintf('Time: %.2f Sec', m/neuron.Fs), 'color', 'w');
    axis equal off tight; title('A\cdot C'); colorbar; 
    
    drawnow();
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [600, 1248]);
        avi_file.writeVideo(temp);
    end
end
if save_avi
    avi_file.close();
end


%% play videos to verify the demixing results  
cell_id = [156, 182, 208, 195, 130, 167, 181];
nc = ceil((length(cell_id)+1)/2);
ctr = round( neuron.estCenter()); 
center = ctr(cell_id, :);
r0 = max(1, min(center(:, 1))-20); r1 = min(d1, max(center(:, 1))+20);
c0 = max(1, min(center(:, 2))-20); c1 = min(d2, max(center(:, 2))+20);
indr = r0:r1;
indc = c0:c1;
center = bsxfun(@minus, center, [r0, c0]-1);
Ysignal = neuron.reshape(Ysignal, 2); 
Ysignal_box = Ysignal(r0:r1, c0:c1, :);

figure('position', [100, 500, 1400, 480]);
avi_file = VideoWriter([dir_nm, file_nm, '_patch_neurons.avi']);
avi_file.open();
ACmax = max(reshape(neuron.A(:, cell_id)*neuron.C(cell_id, :), 1, []))*0.7;
ACmin = 1;
subplot(2, round(nc/2)+1, 1);
h_img = imagesc(Ysignal_box(:, :, 1), [ACmin, ACmax]); hold on;
axis equal off tight;
for m=1:length(cell_id)
    temp = text(center(m,2), center(m, 1), num2str(m), 'color', 'g', 'fontsize', 12, 'fontweight', 'bold');
    tmp_ctr = neuron.Coor{cell_id(m)};
    temp = plot(tmp_ctr(1, 2:end)-c0, tmp_ctr(2, 2:end)-r0, 'r');
end
for t=1:2:T
    subplot(2, nc, 1); hold off;
    imagesc(Ysignal_box(:, :, t), [ACmin, ACmax]);hold on;
    axis equal off tight;
    for m=1:length(cell_id)
        text(center(m,2), center(m, 1), num2str(m), 'color', 'g', 'fontsize', 12, 'fontweight', 'bold');
        tmp_ctr = neuron.Coor{cell_id(m)};
        plot(tmp_ctr(1, 5:end)-c0, tmp_ctr(2, 5:end)-r0, 'r');
    end
    title(sprintf('Time: %.2f Sec', t/neuron.Fs));
    
    for m=1:length(cell_id)
        subplot(2, nc, m+1);
        img = neuron.reshape(neuron.A(:, cell_id(m))*neuron.C(cell_id(m), t), 2);
        img = img(indr, indc);
        imagesc(img, [ACmin, ACmax]);
        axis equal off tight;
        title(sprintf('neuron %d', m));
    end
    drawnow();
    frame = getframe(gcf);
    frame.cdata = imresize(frame.cdata, [480, 1400]);
    avi_file.writeVideo(frame);
    % disp(t);
end
avi_file.close();