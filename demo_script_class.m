clear;
gcp;
% same demo as demo_script.m but using the class @CNMF
%% load file

addpath(genpath('utilities'));
addpath(genpath('deconvolution'));
CNM = CNMF;
CNM.file = 'demoMovie.tif';  % insert path to tiff stack here  

CNM.optionsSet(...                                    % dimensions of datasets
    'search_method','dilate',...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',0.8...                    % merging threshold
    ); 

sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)

CNM.read_file(sframe,num2read)

%todo integrate motion correction
                                     
%% Data pre-processing
CNM.p = 2;
CNM.preprocess;

%% fast initialization of spatial components using greedyROI and HALS
K = 40; tau = 4;  % number of components to be found     % std of gaussian kernel (size of neuron) 
CNM.initComponents(K, tau);

% display centers of found components
% figure;imagesc(Source.Cn);
%     axis equal; axis tight; hold all;
%     scatter(Source.cm(:,2),Source.cm(:,1),'mo');
%     title('Center of ROIs found from initialization algorithm');
%     drawnow;
    
%% update spatial components
CNM.updateSpatial();
%% update temporal components
CNM.updateTemporal(0);

%% classify +traces
CNM.evaluateComponents();
CNM.CNNClassifier('')
%%
CNM.keepComponents()

%% merge found components
Apr = CNM.A;    % store non-merged components
Cpr = CNM.C;
CNM.merge();

display_merging = true; % flag for displaying merging example
try
    i = 1; randi(length(Source.ROI));
catch
    disp('no merged ROIS');
    display_merging = false;
end
if display_merging
    ln = length(Source.ROI{i});
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(Apr(:,CNM.merged_ROIs{i}(j)),CNM.options.d1,CNM.options.d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(CNM.A(:,CNM.K-length(CNM.merged_ROIs)+i),CNM.dims));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:CNM.T,(diag(max(Cpr(CNM.merged_ROIs{i},:),[],2))\Cpr(CNM.merged_ROIs{i},:))'); 
            hold all; plot(1:CNM.T,CNM.C(CNM.K-length(CNM.ROI)+i,:)/max(CNM.C(CNM.K-length(CNM.merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
CNM.updateSpatial();
CNM.updateTemporal(2);
%% extract DF/F
CNM.extractDFF(); % extract DF/F values.

%% do some plotting
%contour_threshold = 0.95;   % amount of energy used for each component to construct contour plot
%figure;
%CNM.viewContours(contour_threshold, 1);

%%
CNM.plotComponentsGUI();     % display all components
%pause;
%% make movie
%Source.makePatchVideo() 