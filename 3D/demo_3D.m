clear;
%% load file (courtesy of D. Pacheco and M. Murthy, Princeton University)

addpath(genpath('../../ca_source_extraction'));
nam = 'data3D.mat';
if ~exist(nam,'file')  % download file if it doesn't exist in the directory
    url = 'https://www.dropbox.com/s/ii3aji4my1i81n2/data3D.mat?dl=1';
    outfilename = websave(nam,url);    
end

load(nam);

if ndims(Y) == 4
    [d1,d2,d3,T] = size(Y);                            % dimensions of dataset
else
    [d1,d2,T] = size(Y);
    d3 = 1;
end
d = d1*d2*d3;                                          % total number of pixels

%% Set parameters

K = 200;                                          % number of components to be found
tau = [3,3,2];                                    % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics for slow imaging rate)
merge_thr = 0.95;                                 % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'d3',d3,...                  % dimensions of datasets
    'search_method','dilate',...                 % search locations when updating spatial components
    'maxIter',15,...                             % number of NMF iterations during initialization
    'deconv_method','constrained_foopsi',...     % activity deconvolution method
    'temporal_iter',2,...                        % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau,'nb',1 ...
    );
%% Data pre-processing

[P,Y] = preprocess_data(Y,p);
Cn = correlation_image_3D(Y); % for large datasets change with reshape(P.sn,d1,d2,d3), %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)

%% show movie of projections
% plot4Dproj(Y, Cn, [d1,d2,d3]);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize
ff = find(sum(Ain)<1e-3*mean(sum(Ain)));   % remove very small components
Ain(:,ff) = [];
Cin(ff,:) = [];
center(ff,:) = [];

%% display centers of found components
plotCenteroverY(Cn, center, [d1,d2,d3]);  % plot found centers against max-projections of background image

%% update spatial components
Yr = reshape(Y,d,T);
%clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%%
%[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Yr,A,C,b,f,YrA,options);
%% plot components
%plot_components_3D_GUI(Y,A,C,b,f,Cn,options)
%% merge found components and repeat (optional)
%[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

%% repeat (optional)
%[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,Pm,options);
%[C2,f2,P2,S2] = update_temporal_components(Yr,A2,b2,Cm,f,Pm,options);

plot_components_3D_GUI(Y,A,C,b,f,Cn,options);