classdef CNMF < handle
   
    % Class for the CNMF object for standard 2p motion correction
    %% properties
    properties
        file = '';              % path to file
        Y;                      % raw data in X x Y ( x Z) x T format
        data;                   % pointer to memory mapped file        
        Yr;                     % raw data in 2d matrix XY(Z) X T format
        sn;                     % noise level for each pixel
        A;                      % spatial components of neurons
        b;                      % spatial components of background
        C;                      % temporal components of neurons
        f;                      % temporal components of backgrounds
        S;                      % neural activity
        Y_res;                  % residual signal (Y - AC - bf)
        R;                      % matrix of component residuals
        bl;                     % baseline of each component
        c1;                     % initial value of each component
        g;                      % discrete time constants for each component
        neuron_sn;              % noise level for each temporal component
        T;                      % number of timesteps
        dims;                   % dimensions of FOV
        d;                      % total number of pixels/voxels
        cm;                     % center of mass for each component
        Coor;                   % neuron contours
        Df;                     % background for each component to normalize the filtered raw data
        C_df;                   % temporal components of neurons and background normalized
        F0;                     % baseline fluorescence for each trace (constant or varying)
        F;                      % raw fluorescence traces (C + R)
        options;                % options for model fitting
        gSig = [5,5];           % half size of neuron
        P;                      % some estimated parameters
        fr = 30;                % frame rate in Hz
        decay_time = 0.4;       % indicator decay time in seconds
        K;                      % number of components
        nb = 1;                 % number of background components
        gnb = 2;                % number of global background components (used only for processing in patches)
        p = 1;                  % order of AR system
        minY;                   % minimum data value
        nd;                     % FOV dimensionality (2d or 3d)
        keep_cnn;               % CNN classifier acceptance
        val_cnn;                % CNN classifier value
        rval_space;             % correlation in space
        rval_time;              % correlation in time
        max_pr;                 % max_probability
        sizeA;                  % size of component
        keep_eval;              % evaluate components acceptance
        fitness;                % exceptionality measure
        keep_exc;               % event exceptionality acceptance 
        thr_exc;                % event exceptionality threshold
        A_throw;                % rejected spatial components
        C_throw;                % rejected temporal components
        S_throw;                % neural activity of rejected traces
        R_throw;                % residuals of rejected components
        A_pre;                  % spatial components before merging
        C_pre;                  % temporal components before merging
        ind_keep;               % selected components (binary vector)
        ind_throw;              % rejected components (binary vector)
        merged_ROIs;            % indeces of merged components
        CI;                     % correlation image
        contours;               % contour coordinates for each spatial component
        patches;                % coordinates of each patch
        patch_size = [64,64];   % size of each patch
        overlap = [10,10];      % overlap between patches
        memmaped = false;       % 
        RESULTS                 % object array with results on each patch
        indicator = 'GCaMP6f';  
        kernel;        
        kernels;
    end
    
    methods
        
        %% construct object and set options
        function obj = CNMF(varargin)
            obj.options = CNMFSetParms();
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:});
            end
        end        
        
        %% read file possibly memory map it and set dimensionality variables
        function readFile(obj,filename,memmaped,sframe,num2read)
            if exist('memmaped','var'); obj.memmaped = memmaped; end
            if ~exist('sframe','var'); sframe = []; end
            if ~exist('num2read','var'); num2read = []; end
            if ischar(filename)
                if obj.memmaped
                    [filepath,name,ext] = fileparts(filename);
                    obj.file = fullfile(filepath,[name,'.mat']);
                    if exist(obj.file,'file')                    
                        obj.Y = matfile(obj.file,'Writable',true);
                    else
                        obj.Y = memmap_file(filename,sframe,num2read,5000);
                    end
                else
                    obj.file = filename;
                    obj.Y = single(read_file(filename,sframe,num2read));
                end                
            else        % filename is an array already loaded in memory
                obj.memmaped = false;
                %obj.Y = single(read_file(filename,sframe,num2read));
                obj.Y = single(filename);            
            end
            
            if obj.memmaped
                obj.minY = obj.Y.nY;
                dimsY = size(obj.Y,'Y');
            else
                obj.minY = min(obj.Y(:));
                obj.Y = obj.Y - obj.minY;
                dimsY = size(obj.Y);
            end                        
            obj.nd = ndims(dimsY)-1;                              
            obj.T = dimsY(end);
            obj.dims = dimsY(1:end-1);
            obj.d = prod(obj.dims);
            obj.options.d1 = dimsY(1);
            obj.options.d2 = dimsY(2);
            if obj.nd > 2; obj.options.d3 = dimsY(3); end
            if obj.memmaped
                obj.Yr = obj.Y;
            else
                obj.Yr = reshape(obj.Y,obj.d,obj.T);
            end
        end
        
        %% create memory mapped file for patch processing
        
        function CNMF_patches(obj,filename,memmaped)
            if exist('memmaped','var'); obj.memmaped = memmaped; end
            if obj.memmaped
                [filepath,name,ext] = fileparts(filename);
                obj.file = fullfile(filepath,[name,'.mat']);
                if exist(obj.file,'file')                    
                    obj.data = matfile(obj.file,'Writable',true);
                else
                    obj.data = memmap_file(filename,1,[],5000);
                end
            else
                obj.file = 'filename';
                obj.data = read_file(filename);
            end
            obj.options = CNMFSetParms();
            sizY = size(obj.data,'Y');
            obj.dims = sizY(1:end-1);
            obj.nd = length(obj.dims);
            obj.minY = obj.data.nY;
            obj.T = sizY(end);
        end
        
        %% create object with pre-loaded array
        function loadArray(obj,Y)
            obj.Y = single(Y);
            obj.minY = min(obj.Y(:));
            obj.Y = obj.Y - obj.minY;                             % make data non-negative
            obj.nd = ndims(obj.Y)-1;                              
            dimsY = size(obj.Y);
            obj.T = dimsY(end);
            obj.dims = dimsY(1:end-1);
            obj.d = prod(obj.dims);
            obj.Yr = reshape(obj.Y,obj.d,obj.T);
            obj.options.d1 = dimsY(1);
            obj.options.d2 = dimsY(2);
            if obj.nd > 2; obj.options.d3 = dimsY(3); end
        end
        
        %% update options
        function optionsSet(obj, varargin)
            obj.options = CNMFSetParms(obj.options, varargin{:});
            obj.p = obj.options.p;
        end
        
        
        %% data preprocessing
        function preprocess(obj)
            [obj.P,obj.Y] = preprocess_data(obj.Y,obj.p,obj.options);
            obj.sn = obj.P.sn;
        end        
        
        %% fast initialization
        function initComponents(obj, K, gSig)
            if exist('gSig','var'); obj.gSig = gSig; end
            [obj.A, obj.C, obj.b, obj.f, obj.cm] = initialize_components(obj.Y, K, obj.gSig, obj.options);
            obj.K = K;
        end
        
        %% update spatial components
        function updateSpatial(obj)        
            if strcmpi(obj.options.spatial_method,'regularized')
                A_ = [obj.A, obj.b];
            else
                A_ = obj.A;
            end
            [obj.A, obj.b, obj.C] = update_spatial_components(obj.Yr, obj.C, obj.f, A_, obj.P, obj.options);
        end
        
        %% update temporal components
        function updateTemporal(obj,p)
            if ~exist('p','var'); p = obj.p; end
            obj.P.p = p;
            [obj.C, obj.f, obj.P, obj.S,obj.R] = update_temporal_components(...
                obj.Yr, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
            obj.bl = cell2mat(obj.P.b);
            obj.c1 = cell2mat(obj.P.c1);
            obj.g = obj.P.gn;
            obj.neuron_sn = cell2mat(obj.P.neuron_sn);
        end
        
        %% deconvolve
        function deconvolve(obj) % TODO: parallel implementation
            if isempty(obj.F)
                if isempty(obj.R)
                     compute_residuals(obj);
                end
                obj.F = obj.C + obj.R;
            end
            if obj.p == 1; model_ar = 'ar1'; elseif obj.p == 2; model_ar = 'ar2'; else; error('non supported AR order'); end
            for i = 1:size(obj.F,1)
                spkmin = GetSn(obj.F(i,:));
                [cc, spk, opts_oasis] = deconvolveCa(obj.F(i,:),model_ar,'optimize_b',true,'method','thresholded',...
                    'optimize_pars',true,'maxIter',20,'smin',spkmin,'window',200); 
                obj.bl = opts_oasis.b;
                obj.C(i,:) = full(cc(:)' + cb);
                obj.S(i,:) = spk(:)';
                obj.R(i,:) = obj.F(i,:) - obj.C(i,:);
                obj.c1(i) = 0;
                obj.neuron_sn(i) = opts_oasis.sn;
                obj.g{i} = opts_oasis.pars(:)';
            end
        end
        
        %% merge components
        function merge(obj)
            obj.A_pre = obj.A;
            obj.C_pre = obj.C;
            [obj.A, obj.C, obj.K, obj.merged_ROIs, obj.P, obj.S] = merge_components(...
                obj.Yr,obj.A, [], obj.C, [], obj.P,obj.S, obj.options);
        end 
        
        %% extract DF_F
        function extractDFF(obj)
            %[obj.C_df,~] = extract_DF_F(obj.Yr,obj.A,obj.C,obj.P,obj.options);
            [obj.C_df,obj.F0] = detrend_df_f(obj.A,obj.b,obj.C,obj.f,obj.R,obj.options);
        end

        %% CNN classifier
        function CNNClassifier(obj,classifier)
            try
                [obj.keep_cnn,obj.val_cnn] = cnn_classifier(obj.A,obj.dims,classifier,obj.options.cnn_thr);
            catch
                obj.keep_cnn = true(size(obj.A,2),1);
                obj.val_cnn = ones(size(obj.A,2),1);
            end
        end
        
        %% evaluate components
        function evaluateComponents(obj)
            [obj.rval_space,obj.rval_time,obj.max_pr,obj.sizeA,obj.keep_eval] = classify_components(...
                obj.Y,obj.A,obj.C,obj.b,obj.f,obj.R,obj.options);
        end

        %% keep components
        function keepComponents(obj,ind_keep)
            if ~exist('ind_keep','var')
                obj.ind_keep = (obj.keep_eval | obj.keep_cnn) & obj.keep_exc;
            else
                obj.ind_keep = ind_keep;
            end
            obj.ind_throw = ~obj.ind_keep;
            obj.A_throw = obj.A(obj.ind_throw,:);
            obj.C_throw = obj.C(obj.ind_throw,:);
            obj.S_throw = obj.S(obj.ind_throw,:);
            obj.R_throw = obj.R(obj.ind_throw,:);
            obj.A = obj.A(:,obj.ind_keep);
            obj.C = obj.C(obj.ind_keep,:);
            obj.S = obj.S(obj.ind_keep,:);
            obj.R = obj.R(obj.ind_keep,:);
            obj.bl = obj.bl(obj.ind_keep);
            obj.c1 = obj.c1(obj.ind_keep);
            obj.neuron_sn = obj.neuron_sn(obj.ind_keep);
            obj.g = obj.g(obj.ind_keep);
        end
        
        %% normalize components
        function normalize(obj)
            nA = sqrt(sum(obj.A.^2,1));
            obj.A = bsxfun(@times, obj.A, 1./nA);
            obj.C = bsxfun(@times, obj.C, nA');
            if ~isempty(obj.S)
                obj.S = bsxfun(@times, obj.S, nA');
            end
            if ~isempty(obj.R)
                obj.R = bsxfun(@times, obj.R, nA');
            end
            nB = sqrt(sum(obj.b.^2,1));
            obj.b = bsxfun(@times, obj.b, 1./nB);
            obj.f = bsxfun(@times, obj.f, nB');
        end
           
        %% compute residuals
        function compute_residuals(obj)
            AA = obj.A'*obj.A;
            AY = mm_fun(obj.A,obj.Y);
            nA2 = sum(obj.A.^2,1);
            obj.R = bsxfun(@times, AY - AA*obj.C - (obj.A'*obj.b)*obj.f,1./nA2);
        end
        
        %% compute event exceptionality
        function eventExceptionality(obj)
            if isempty(obj.R)
                compute_residuals(obj);
            end
            obj.fitness = compute_event_exceptionality(obj.C + obj.R,obj.options.N_samples_exc,obj.options.robust_std);
            obj.keep_exc = (obj.fitness < obj.options.min_fitness);
        end
        
        %% correlation image
        function correlationImage(obj,sz,batch_size,min_batch_size)
            if ~exist('sz','var'); sz = []; end
            if ~exist('batch_size','var'); batch_size = []; end
            if ~exist('min_batch_size','var'); min_batch_size = []; end
            obj.CI = correlation_image_max(obj.Y,sz,obj.dims,batch_size,min_batch_size);
        end
        
        %% center of mass
        function COM(obj)
            if obj.nd == 2
                obj.cm = com(obj.A,obj.dims(1),obj.dims(2));
            elseif obj.nd == 3
                obj.cm = com(obj.A,obj.dims(1),obj.dims(2),obj.dims(3));
            end                
        end
        
        %% plot contours
        function plotContours(obj,display_numbers,max_number,ln_wd,ind_show)
            if ~exist('display_numbers','var'); display_numbers = []; end
            if ~exist('max_number','var'); max_number = []; end
            if ~exist('ln_wd','var'); ln_wd = []; end
            if ~exist('ind_show','var'); ind_show = []; end
            if isempty(obj.CI); Cn = reshape(obj.sn,obj.dims); else; Cn = obj.CI; end
            [obj.contours,~,~] = plot_contours(obj.A,Cn,obj.options,display_numbers,max_number,obj.contours, ln_wd, ind_show,obj.cm);
        end
        
        %% plot components GUI
        function plotComponentsGUI(obj)
            if or(isempty(obj.CI), ~exist('CI', 'var') )
                correlationImage(obj);
            end
            plot_components_GUI(obj.Yr,obj.A,obj.C,obj.b,obj.f,obj.CI,obj.options)
        end
        
        %% plot centers
        function plotCenters(obj)
            if isempty(obj.CI); Cn = reshape(obj.sn,obj.dims); else; Cn = obj.CI; end
            COM(obj);
            figure;imagesc(Cn); axis equal; axis tight; hold all;
                scatter(obj.cm(:,2),obj.cm(:,1),'mo');
                title('Center of ROIs found from initialization algorithm');
                drawnow; hold off;
        end
        
        %% display merging
        function displayMerging(obj)
            display_merging = true; % flag for displaying merging example
            try
                i = 1; randi(length(obj.merged_ROIs));
            catch
                disp('no merged ROIS');
                display_merging = false;
            end
            if display_merging
                ln = length(obj.merged_ROIs{i});
                figure;
                    set(gcf,'Position',[300,300,(ln+2)*300,300]);
                    for j = 1:ln
                        subplot(1,ln+2,j); imagesc(reshape(obj.A_pre(:,obj.merged_ROIs{i}(j)),obj.options.d1,obj.options.d2)); 
                            title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
                    end
                    subplot(1,ln+2,ln+1); imagesc(reshape(obj.A(:,obj.K-length(obj.merged_ROIs)+i),obj.dims));
                            title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
                    subplot(1,ln+2,ln+2);
                        plot(1:obj.T,(diag(max(obj.C_pre(obj.merged_ROIs{i},:),[],2))\obj.C_pre(obj.merged_ROIs{i},:))'); 
                        hold all; plot(1:obj.T,obj.C(obj.K-length(obj.merged_ROIs)+i,:)/max(obj.C(obj.K-length(obj.merged_ROIs)+i,:)),'--k')
                        title('Temporal Components','fontsize',16,'fontweight','bold')
                    drawnow;
            end
        end
        
        %% fit transform
        function fit(obj,Y,options,K,gSig)
            if exist('gSig','var'); obj.gSig = gSig; end
            obj.readFile(Y);
            obj.optionsSet(options);       % set up options parameters
            obj.preprocess();
            obj.initComponents(K);
            obj.updateSpatial();
            obj.updateTemporal(0);
            obj.evaluateComponents();
            obj.CNNClassifier('cnn_model.h5');
            obj.eventExceptionality();
            obj.keepComponents();
            obj.merge();
            obj.updateSpatial();
            obj.updateTemporal(obj.p);
        end
        
        
        %% construct patches
        function createPatches(obj,patch_size,overlap)
            if exist('patch_size','var'); obj.patch_size = patch_size; end
            if exist('overlap','var');  obj.overlap = overlap; end
            obj.patches = construct_patches(obj.dims,obj.patch_size,obj.overlap);
        end
        
        %% fit patches
        function fitPatches(obj)
            obj.options.refine_flag = false;
            [obj.A,obj.b,obj.C,obj.f,obj.S,Pr,obj.RESULTS,obj.R] = ... 
                run_CNMF_patches(obj.Y,obj.K,obj.patches,obj.gSig,obj.p,obj.options);
            obj.bl = cell2mat(Pr.b);
            obj.c1 = cell2mat(Pr.c1);
            obj.neuron_sn = cell2mat(Pr.neuron_sn);
            obj.g = cell2mat(Pr.gn);
            obj.sn = Pr.sn;
            obj.P.sn = Pr.sn;
        end
        
    end
end