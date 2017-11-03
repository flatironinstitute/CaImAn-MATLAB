clear;
%% load file
addpath(genpath('../ca_source_extraction'));  

clear;
gcp;    % start a local cluster

filename = 'demoSue2x.tif';
if ~exist(filename,'file');
    url = 'https://www.dropbox.com/s/36xdfd28eone0hj/demoSue2x.tif?dl=1';
    fprintf('downloading the file...');
    outfilename = websave(filename,url);
    fprintf('done. \n');
end


Y = read_file(filename);
Y = Y - min(Y(:)); 
if ~isa(Y,'single');  Y = single(Y);  end           % convert to single
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% Set parameters
% All the parameters can be changed here look at the CNMFsetparams.m file 
% for more information about each parameters.
K = 130;                                           % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','dilate','dist',3,...       % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );
%% Data pre-processing

[P,Y] = preprocess_data(Y,p);
[A,C,b,f,~] = initialize_components(Y,K,tau,options,P);  % initialize
Cn =  correlation_image(Y);
%% GUI ------
%this will create a User interface for the rest of the pipeline described in demoscript and 
%demo pipeline. LOOK at demo pipeline to understand how to use this GUI and what to pass it
%           
%It is a mix of everyhting found in those pipelines. In the future it will be also called by 
%the initialization GUI 
%
%run GUI will let you, refine the components manually, change the parameters and see the results,
%analyse the trace of each components, use a classification algorithm to
%find the components instead, add and remove components, save the ROIS
%and finish the pipeline with button clicks. 
%It is still an optional method.

% WARNING : do not use simplemode for now - minimum MATLAB VERSION 2014
%                   save is only compatible with the JSON matlabpackage

% here is what the GUI needs to receive in parameters
GUIout = ROI_GUI(Y,A,P,options,Cn,C,b,f);   
pause;
%% -----------
A =GUIout{1}; 
options = GUIout{2};
Cdec = GUIout{3};
ROIvars = GUIout{4};
YrA=GUIout{5};
Contours=GUIout{6};
b=GUIout{7};
f=GUIout{8};
Spikes=GUIout{9};

%% make movie

make_patch_video(A,ROIvars.C,b,f,YrA,Contours,options)