clear;
% same demo as demo_script.m but using the class @Sources2D
% (possibly outdated at the moment)
%% load file

addpath(genpath('utilities'));
addpath(genpath('deconvolution'));
Source = Sources2D;
Source.file = 'demoMovie.tif';  % insert path to tiff stack here  

Source.updateParams(...                                    % dimensions of datasets
    'search_method','dilate',...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',0.8...                    % merging threshold
    ); 

sframe=1;						% user input: first frame to read (optional, default 1)
num2read=2000;					% user input: how many frames to read   (optional, default until the end)
Source.read(sframe,num2read)

%%%% MOTION CORRECTION
%todo integrate motion correction
                                     
%% Data pre-processing
p = 2;    % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and deca 
Source.preprocess(p);

%% fast initialization of spatial components using greedyROI and HALS
K = 40; tau = 4;  % number of components to be found     % std of gaussian kernel (size of neuron) 
Source.initComponents(K, tau);

% display centers of found components
figure;imagesc(Source.Cn);
    axis equal; axis tight; hold all;
    scatter(Source.cm(:,2),Source.cm(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%%%%%%%%%%%%%%%introduction to GUI mode ;
%% run GUI will let you, refine the components manually, change the parameters and see the results,
% analyse the trace of each components, use a classification algorithm to
% find the components instead, add and remove components, save the ROIS
% and finish the pipeline with button clicks. 

%it is still an optional method.
run_GUI = false;
if run_GUI
    %todo to implement 
    ROIvars.C = obj.C;
    Coor = plot_contours(obj.A,obj.Cn,obj.options,1); close;
    % here is what the GUI needs to receive in parameters
    GUIout = ROI_GUI(obj.Y ,obj.A ,obj.P ,obj.options ,obj.Cn ,Coor ,...
        obj.keep ,ROIvars ,obj.b ,obj.f ,obj.S,obj.Yra);   
    options = GUIout{2};
    keep = GUIout{3};    
end

%%OR
%% manually refine components (optional)
refine_components = false;  % flag for manual refinement
if refine_components
    Source.refineComponents();
end
    
%% update spatial components
Source.updateSpatial();
%% update temporal components
Source.updateTemporal();

%%%%%%%%%%%%%%%% classify +traces
Source.classifyComponents();
Source.compute_event_exceptionality();
%% merge found components
Apr = Source.A;    % store non-merged components
Cpr = Source.C;
Source.merge();

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
            subplot(1,ln+2,j); imagesc(reshape(Apr(:,Source.ROI{i}(j)),Source.options.d1,Source.options.d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Source.A(:,Source.K-length(Source.ROI)+i),Source.options.d1,Source.options.d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:Source.T,(diag(max(Cpr(Source.ROI{i},:),[],2))\Cpr(Source.ROI{i},:))'); 
            hold all; plot(1:Source.T,Source.C(Source.K-length(Source.ROI)+i,:)/max(Source.C(Source.K-length(Source.ROI)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
Source.updateSpatial();
Source.updateTemporal();
%% do some plotting
Source.orderROIs();     % order components
Source.extractDF_F(); % extract DF/F values.

%%%%%%%%%%%%%% deconvolve
%Cdec = Source.deconvTemporal();

contour_threshold = 0.95;   % amount of energy used for each component to construct contour plot
figure;
Source.viewContours(contour_threshold, 1);

%%%%%%% save
save = false;
if save
    savejson('jmesh',Source.json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
end
Source.plotComponentsGUI();     % display all components
pause;
%% make movie
Source.makePatchVideo() 