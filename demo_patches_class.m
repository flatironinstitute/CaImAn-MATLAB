% demo script for splitting the field of view in patches and processing in 
% parallel with or without memory mapping. The demo follows the file 
% demo_patches.m using the CNMF class object

clear;
%% setup path to file and package
gcp;
path_to_package = '../ca_source_extraction';   % path to the folder that contains the package
addpath(genpath(path_to_package));
             
filename = '/Users/epnevmatikakis/Documents/Ca_datasets/Neurofinder/neurofinder.02.00/images/neurofinder0200_rig.tif';      
        % path to file (assumed motion corrected)
        
is_memmaped = true;        % choose whether you want to load the file in memory or not

%% create object and load file

CNM = CNMF();
if is_memmaped
    CNM.readFile(filename,is_memmaped);
else
    CNM.readFile(filename,is_memmaped,1,2000); % load only a part of the file due to memory
end

%% set options and create patches

patch_size = [32,32];    % size of each patch along each dimension (optional, default: [32,32])
overlap = [6,6];         % amount of overlap in each dimension (optional, default: [4,4])
K = 10;                  % number of components to be found
gSig = 7;                % std of gaussian kernel (size of neuron) 
p = 2;                   % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
gnb = 3;                 % order of background
merge_thr = 0.8;         % merging threshold

options = CNMFSetParms(...
    'd1',CNM.dims(1),'d2',CNM.dims(2),...
    'search_method','dilate',...                % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'nb',1,...                                  % number of background components per patch
    'gnb',gnb,...                               % number of global background components
    'ssub',2,...
    'tsub',1,...
    'p',p,...                                   % order of AR dynamics
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',gSig,... 
    'spatial_method','regularized',...
    'cnn_thr',0.2,...
    'patch_space_thresh',0.25,...
    'min_SNR',2);

CNM.optionsSet(options);
CNM.gnb = gnb;
CNM.K = K;
CNM.patch_size = patch_size;                % size of each patch along each dimension (optional, default: [32,32])
CNM.overlap = overlap;                      % amount of overlap in each dimension (optional, default: [4,4])
CNM.createPatches();                        % create patches


%% fit all patches
CNM.fitPatches();

%% component classification

CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
CNM.CNNClassifier('')       % evaluate spatial components with the CNN classifier
CNM.eventExceptionality();  % evaluate traces
CNM.keepComponents();       % keep the components that are above certain thresholds

%% repeat processing

CNM.updateSpatial();
CNM.updateTemporal();
CNM.extractDFF();            % extract DF/F values.

%% do some plotting
figure;
CNM.correlationImage();
CNM.plotContours();
CNM.plotComponentsGUI();     % display all components