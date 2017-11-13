clear;
gcp;

% same demo as demo_script.m but using the class @CNMF
%% load file

addpath(genpath('utilities'));
addpath(genpath('deconvolution'));
CNM = CNMF;                                     % contruct CNMF object
filename = 'demoMovie.tif';                     % filename to be processed
K = 40;                                         % number of components to be found

options = CNMFSetParms(...   
    'p',2,...                                   % order of AR dynamics    
    'gSig',5,...                                % half size of neuron
    'merge_thr',0.80,...                        % merging threshold  
    'nb',2,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    'space_thresh',0.5,...                      % space correlation threshold
    'cnn_thr',0.2...                            % threshold for CNN classifier    
    );

%%
% Below is the standard processing pipeline. This processing can be
% executed in one shot using the CNM.fit function:
% CNM.fit(filename,options,K)


%% load the dataset and create the object
CNM.readFile(filename);                         % insert path to file here  
CNM.optionsSet(options);                        % setup the options structure

%% Process the dataset

CNM.preprocess;             % preprocessing (compute some quantities)
CNM.initComponents(K);      % initialization
CNM.plotCenters()           % plot center of ROIs detected during initialization
CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)

%% component classification

CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
CNM.CNNClassifier('')       % evaluate spatial components with the CNN classifier
CNM.eventExceptionality();  % evaluate traces
CNM.keepComponents();       % keep the components that are above certain thresholds

%% merge found components
CNM.merge();
CNM.displayMerging();

%% repeat processing

CNM.updateSpatial();
CNM.updateTemporal();
CNM.extractDFF();           % extract DF/F values.

%% do some plotting
figure;
CNM.plotContours();
CNM.plotComponentsGUI();     % display all components