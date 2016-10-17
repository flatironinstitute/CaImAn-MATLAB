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
neuron_raw.Fs = 10;         % frame rate
neuron_raw.options.Fs = 20; 
ssub = 1;           % spatial downsampling factor 
tsub = 1;           % temporal downsampling factor 
neuron_raw.updateParams('ssub', ssub,...  % spatial downsampling factor
    'tsub', tsub, ...  %temporal downsampling factor
    'gSig', 4,... %width of the gaussian kernel, which can approximates the average neuron shape
    'gSiz', 15, ...% average size of a neuron
    'bSiz', 2, ... % half width of the kernel to dilate nonzero pixels of each neuron
    'search_method', 'dilate', ... % searching method
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
debug_on = false;        % debug mode 
save_avi = false;
neuron.options.min_corr = 0.8;
neuron.options.nk = 1; %round(T/(60*neuron.Fs)); % number of knots for spline basis, the interval between knots is 180 seconds
patch_par = [3, 3];  % divide the optical field into 3 X 3 patches and do initialization patch by patch
K = 500; % maximum number of neurons to search within each patch. you can use [] to search the number automatically
neuron.options.bd = []; % boundaries to be removed due to motion correction 
[center, Cn, pnr] = neuron.initComponents_endoscope(Y, K, patch_par, debug_on, save_avi); 
neuron_init = neuron.copy(); 
[~, srt] = sort(max(neuron.C, [], 2)./get_noise_fft(neuron.C), 'descend'); 
neuron.orderROIs(srt); 

%% merge neurons, order neurons and delete some low quality neurons (cell 0, before running iterative udpates)
[Ain, Cin] = neuron.snapshot();   % keep the initialization results
neuron.options.merge_thr = 0.7;     % threshold for merging neurons
merged_ROI = neuron.quickMerge('C');  %merge neurons based on their temporal correlation 
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

%% udpate background (cell 1, the following three blocks can be run iteratively)
% determine nonzero pixels for each neuron 
max_min_ratio = 50;     % it thresholds the nonzero pixels to be bigger than max / max_min_ratio.
neuron.trimSpatial(max_min_ratio);
neuron.options.se = strel('disk', 10); % if neuron contains dendrites and axons, use large values (e.g., 20), otherwise <5 is OK
IND = determine_search_location(neuron.A, 'dilate', neuron.options);

% start approximating theb background
tic;
Ybg = Y-neuron.A*neuron.C; 
ssub = 3;   % downsample the data to improve the speed 
rr = neuron.options.gSiz;  % average neuron size, it will determine the neighbors for regressing each pixel's trace
active_px = (sum(IND, 2)>0);  %If some missing neurons are not covered by active_px, use [] to replace IND
Ybg = neuron.localBG(Ybg, ssub, rr, active_px); % estiamte local background.
fprintf('Time cost in estimating the background:        %.2f seconds\n', toc); 

% subtract the background from the raw data. 
Ysignal = Y - Ybg;
if ~isfield(neuron.P, 'sn')
    neuron.preprocess(Ysignal, 2);
end

%% update spatial components (cell 2), we can iteratively run cell 2& 3 for few times and then run cell 1
spatial_method = 'hals'; % methods for updating spatial components {'hals', 'lars'}, 
% hals is fast, there is no sparse constraint; lars is very slow, but the
% results are sparse. LARS is not stable yet. 
tic;
% update spatial components with model Y = A*C
if strcmpi(spatial_method, 'lars')
    if ~isfield(neuron.P, 'sn')
        neuron.preprocess(Ysignal, 2); 
    end
    neuron.updateSpatial_nb(Ysignal); 
else
    IND = determine_search_location(neuron.A, 'dilate', neuron.options);
    neuron.A = HALS_spatial(Ysignal, neuron.A, neuron.C, IND, 20);
%     neuron.post_process_spatial(); % uncomment this line to postprocess
%     the results 
end
fprintf('Time cost in updating neuronal spatial components:     %.2f seconds\n', toc);

%% update C  (cell 3)
temporal_method = 'hals'; % methods for updating temporal components {'hals', 'deconv'}
tic;
C0 = neuron.C;
if strcmpi(temporal_method, 'deconv')
    neuron.updateTemporal_nb(Ysignal); 
else
    neuron.C = HALS_temporal(Ysignal, neuron.A, neuron.C, 10);
end

fprintf('Time cost in updating neuronal temporal components:     %.2f seconds\n', toc);

%% pick neurons from the residual (cell 3.5). It's not alway necessary 
neuron.manual_add(Ysignal-neuron.A*neuron.C); 

%% display neurons
figure;
neuron.viewNeurons();

%% display contours of the neurons 
figure; 
neuron.viewContours(Cn, 0.9, 0); 
colormap winter; 
axis equal; axis off;
title('contours of estimated neurons');

%% check results by visualizing all spatial and temporal components one by one
folder_nm = [];%'neurons';
neuron.viewNeurons([], C, folder_nm);

%% check spatial and temporal components by playing movies
save_avi = true;
avi_name = 'play_movie.avi';
neuron.runMovie(Ysignal, [0, 100], save_avi, avi_name);
%%
