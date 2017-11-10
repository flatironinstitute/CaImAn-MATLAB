classdef Sources2D < handle
    
    % This class is a wrapper of Constrained NMF for standard 2D data.
    % Author: Pengcheng Zhou, zhoupc1988@gmail.com with modifications from
    % Eftychios Pnevmatikakis
    
    %% properties
    properties
        Y;
        A;          % spatial components of neurons
        T;         
        C;          % temporal components of neurons
        b;          % spatial components of backgrounds
        f;          % temporal components of backgrounds
        S;          % spike counts
        Yr;
        YrA;
        d;
        C_raw;      % raw traces of temporal components
        Cn;         % correlation image
        cm;
        Coor;       % neuron contours
        Df;         % background for each component to normalize the filtered raw data
        C_df;       % temporal components of neurons and background normalized bobj.y Df
        S_df;       % spike counts of neurons normalized bobj.y Df
        options;    % options for model fitting
        tau;
        P;          % some estimated parameters
        Fs = nan;    % frame rate
        indicator = 'GCaMP6f';
        kernel;
        file = '';
        ROI;
        K;
        json_file;
        Ybg;
        Y_res;
        rval_space;
        rval_time;
        max_pr;
        sizeA;
        keep;
        fitness;
        traces;
        fitness_delta;
        Cdec; 
        Pdec; 
        Sdec;
        srt;
        kernels;
    end


    
    %% methods
    methods
        %% constructor and options setting
        function obj = Sources2D(varargin)
            obj.options = CNMFSetParms();
            obj.P = struct('p', 2);
            if nargin>0
                obj.options = CNMFSetParms(obj.options, varargin{:});
            end
            %obj.kernel = create_kernel('exp2');
        end
        
        function read(obj,sframe,num2read)
            obj.Y = bigread2(obj.file,sframe,num2read);
            if ~isa(obj.Y,'double');    
                obj.Y = double(obj.Y);  
            end        
            obj.Y = obj.Y - min(obj.Y(:));        % make data non-negative
            [obj.options.d1,obj.options.d2,obj.T] = size(obj.Y);  % dimensions of dataset
            obj.d = obj.options.d1*obj.options.d2;
            obj.Yr = reshape(obj.Y,obj.d,obj.T);
        end
        

        %% update parameters
        function updateParams(obj, varargin)
            obj.options = CNMFSetParms(obj.options, varargin{:});
        end
        

        %% data preprocessing
        function preprocess(obj,p)
            [obj.P,obj.Y] = preprocess_data(obj.Y,p,obj.options);
        end
        
        
        %% fast initialization
        function initComponents(obj, K, tau)
            if nargin<4 ;    tau = [];             end
            [obj.A, obj.C, obj.b, obj.f, obj.cm] = initialize_components(obj.Y, K, tau, obj.options);
            obj.Cn =  reshape(obj.P.sn,obj.options.d1,obj.options.d2); 
            %correlation_image(obj.Y); %max(obj.Y,[],3); %std(obj.Y,[],3); %
            % image statistic (onlobj.y for displaobj.y purposes)
            obj.tau;
            obj.K;
        end
        

        %% fast initialization for microendoscopic data
                %%TODO EDIT
        [center, Cn, pnr] = initComponents_endoscope(obj, Y, K, patch_sz, debug_on, save_avi);
        

        %% update spatial components
        function updateSpatial(obj)
        % add b to A_, otherwise error in update_spatial_components line 73
        if strcmpi(obj.options.spatial_method,'regularized')
            A_ = [obj.A, obj.b];
        else
            A_ = obj.A;
        end
            [obj.A, obj.b, obj.C] = update_spatial_components(obj.Yr, ...
                obj.C, obj.f, A_, obj.P, obj.options);
        end
        
        
        %% udpate spatial components without background
        function updateSpatial_nb(obj)
            [obj.A, obj.C] = update_spatial_components_nb(obj.Yr, ...
                obj.C, obj.A, obj.P, obj.options);
        end
        
        %% udpate spatial components using NNLS and thresholding
        function updateSpatial_nnls(obj)
            [obj.A, obj.C] = update_spatial_nnls(obj.Yr, ...
                obj.C, obj.A, obj.P, obj.options);
        end
        
        
        %% update temporal components for endoscope data
                %% TODO to edit
        updateSpatial_endoscope(obj, Y, numIter, method, IND_thresh);
        
        
        %% update temporal components
        function updateTemporal(obj)
            [obj.C, obj.f, obj.P, obj.S,obj.YrA] = update_temporal_components(...
                obj.Yr, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        
        
        %% refine components
        function center = refineComponents(obj)
            [obj.A,obj.C,obj.cm] = manuallobj.y_refine_components(...
                obj.Yr,obj.A,obj.C,obj.cm,obj.Cn,obj.tau,obj.options);
            center = obj.cm;
        end
        

        %% udpate temporal components with fast deconvolution
                        %% TODO to edit
        updateTemporal_endoscope(obj, Y, smin)
        

        %% update temporal components without background
        function updateTemporal_nb(obj)
            [obj.C, obj.P, obj.S] = update_temporal_components_nb(...
                obj.Yr, obj.A, obj.b, obj.C, obj.f, obj.P, obj.options);
        end
        

        %% merge found components
        function merge(obj)
            [obj.A, obj.C, obj.K, obj.ROI, obj.P, obj.S] = merge_components(...
                obj.Yr,obj.A, [], obj.C, [], obj.P,obj.S, obj.options);
        end
       

        %% compute the residual
        function [Y_res] = residual(obj)
            obj.Y_res = obj.Yr - obj.A*obj.C - obj.b*obj.f;
            Y_res = obj.Y_res;
        end
       

        %% take the snapshot of current results
        function [A, C,  b, f, P, S] = snapshot(obj)
            A = obj.A;
            C = obj.C;
            b = obj.b;
            f = obj.f;
            P = obj.P;
            try
                S = obj.S;
            catch
                S = [];
            end
        end
       

        %% extract DF/F signal after performing NMF
        function [C_df, Df] = extractDF_F(obj)
            [obj.C_df, obj.Df] = extract_DF_F(obj.Yr, obj.A,obj.C, obj.P, obj.options);       
            C_df =  obj.C_df;
            Df = obj.Df;                  
        end
        

        %% order_ROIs
        function orderROIs(obj, srt)
            %% order neurons
            % srt: sorting order
            nA = sqrt(sum(obj.A.^2));
            nr = length(nA);
            if nargin<2; srt=[]; end
            [obj.A, obj.C, obj.S, obj.P, obj.srt] = order_ROIs(obj.A, obj.C,...
                obj.S, obj.P, srt);
            
            if ~isempty(obj.C_raw)
                if isempty(obj.srt)
                    obj.C_raw = spdiags(nA(:),0,nr,nr)*obj.C_raw;
                end
                obj.C_raw = obj.C_raw(obj.srt, :);
            end
        end
        

        %% view contours
        function viewContours(obj,contour_threshold, display, ind, ln_wd)
            if or(isempty(obj.Cn), ~exist('Cn', 'var') )
                obj.Cn = reshape(obj.P.sn, obj.options.d1, obj.options.d2);
            end
            if (~exist('displaobj.y', 'var') || isempty(display)); displaobj.y=0; end
            if (~exist('ind', 'var') || isempty(ind)); ind=1:size(obj.A, 2); end
            if (~exist('ln_wd', 'var') || isempty(ln_wd)); ln_wd = 1; end
            [obj.Coor, obj.json_file] = plot_contours(obj.A(:, ind), obj.Cn, ...
                contour_threshold, displaobj.y, [], [], ln_wd);
            obj.Coor;
        end
        

        %% plot components
        function plotComponents(obj)
            if or(isempty(obj.Cn), ~exist('Cn', 'var') )
                obj.Cn = [];
            end
            view_components(obj.Yr, obj.A, obj.C, obj.b, obj.f, obj.Cn, obj.options);
        end
        

        %% plot components GUI
        function plotComponentsGUI(obj)
            if or(isempty(obj.Cn), ~exist('Cn', 'var') )
                obj.Cn = [];
            end
            plot_components_GUI(obj.Yr,obj.A,obj.C,obj.b,obj.f,obj.Cn,obj.options)
        end
        

        %% make movie
        function makePatchVideo(obj)
            make_patch_video(obj.A, obj.C, obj.b, obj.f, obj.Yr, obj.Coor,...
                obj.options);
        end

        function classifyComponents(obj)
            [obj.rval_space,obj.rval_time,obj.max_pr,obj.sizeA,obj.keep] = classify_components(...
                obj.Y,obj.A,obj.C,obj.b,obj.f,obj.YrA,obj.options);
        end
        

        %% new methods added bobj.y PC, Since 02/05/2016
        %% down sample data and initialize it
        %TODO to modify 
        [obj_ds, Y_ds] = downSample(obj, Y)
        

        %% up-sample the results
        %TODO to modify 
        upSample(obj, obj_ds, T);
        

        function compute_event_exceptionality(obj)
             obj.traces = prctfilt(obj.C+obj.YrA,8,1000,100);
             obj.fitness = compute_event_exceptionality(obj.traces,0);
             obj.fitness_delta = compute_event_exceptionality(diff(obj.traces,[],2),0);  
        end


        %% copobj.y the objects
        function obj_new = copy(obj)
            % Instantiate new object of the same class.
            obj_new = feval(class(obj));
            % Copobj.y all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                obj_new.(p{i}) = obj.(p{i});
            end
        end
        

        %% quick merge neurons based on spatial and temporal correlation
        %TODO to modify 
        [merged_ROIs, newIDs] = quickMerge(obj, merge_thr);
        
        
        %% quick view
        %TODO to modify 
        viewNeurons(obj, ind, C2, folder_nm);
        displaobj.yNeurons(obj, ind, C2, folder_nm);
        
        %% delete neurons
        function delete(obj, ind)
            obj.A(:, ind) = [];
            obj.C(ind, :) = [];
            if ~isempty(obj.S); obj.S(ind, :) = []; end
            if ~isempty(obj.C_raw); obj.C_raw(ind, :) = []; end
            if isfield(obj.P, 'kernel_pars')&&(  ~isempty(obj.P.kernel_pars))
                obj.P.kernel_pars(ind, :) = [];
            end
        end
        
        %% estimate neuron centers
        function estCenter(obj)
            obj.cm = com(obj.A, obj.options.d1, obj.options.d2);
        end
        

        %% update A & C using HALS
        function [obj, IDs] = HALS_AC(obj)
            %update A,C,b,f with HALS
            obj.Y = obj.reshape(obj.Y, 1);
            [obj.A, obj.C, obj.b, obj.f, IDs] = HALS_2d(obj.Y, obj.A, obj.C, obj.b,...
                obj.f, obj.options);
        end
        
        %% update backgrounds
        %[B, F] = updateBG(obj, obj.Y, nb, model)
        
        %% reshape data
        function Y = reshape(obj, dim)
            % reshape the imaging data into diffrent dimensions
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            if dim==1; Y=reshape(obj.Y, d1*d2, []);  %each frame is a vector
            else  Y = reshape(full(obj.Y), d1, d2, []);    %each frame is an image
            end
        end

        
        %% deconvolve all temporal components
        function Cdec = deconvTemporal(obj)
            Cdec = zeros(size(obj.C));
            kernels = cell(obj.K,1);
            for i = 1:obj.K
                [Cdec(i,:),obj.S(i,:),kernels{i}] = deconvCa(...
                    obj.C(i,:), [], 3, true, false, [], 20, [], 0);
            end
            obj.Cdec = Cdec;
            obj.kernels = kernels;
        end

        
        %% update background
        function [ind_bg, bg] = linearBG(obj)
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            % model the background as linear function
            if ndims(obj.Y)==3; obj.Y = neuron.reshape(obj.Y,1); end
            obj.T = size(obj.Y,2);
            % find background area
            ind_frame = round(linspace(1, obj.T, min(obj.T, 500)));
            tmp_C1 = correlation_image(obj.Y(:, ind_frame), [1, 2], d1, d2);
            tmp_C2 = correlation_image(obj.Y(:, ind_frame), [0,1]+ obj.options.gSiz, d1, d2);
            tmp_Cn = tmp_C1(:)-tmp_C2(:);
            ind = (tmp_Cn>0);
            ind_bg = and(ind, tmp_Cn<quantile(tmp_Cn(ind), 0.1));
            % get the mean activitobj.y within the selected area
            obj.Ybg = mean(obj.Y(ind_bg(:), :), 1);
            f1 = (obj.Ybg-mean(obj.Ybg))/std(obj.Ybg);
            
            % regress Dobj.Y over df to remove the trend
            df = diff(f1,1);
            b1 = diff(obj.Y, 1,2)*df'/(df*df');
            
            Yres = obj.Y-b1*f1;
            b0 = median(obj.Yres, 2);
            %             b0 = quantile(obj.Yres, 0.05,2);
            obj.Y = bsxfun(@minus, Yres, b0);
            bg.b = [b1, b0];
            bg.f = [f1; ones(1,obj.T)];
        end

        
        %% select background pixels
        function [f1, ind_bg] = findBG(obj, q)
            d1 = obj.options.d1;
            d2 = obj.options.d2;
            if nargin<3;    q = 0.1; end;   %quantiles for selecting background pixels
            % model the background as linear function
            if ndims(obj.Y)==3; obj.Y = neuron.reshape(obj.Y,1); end
            obj.T = size(obj.Y,2);
            % find background area
            ind_frame = round(linspace(1, obj.T, min(obj.T, 500)));
            tmp_C1 = correlation_image(obj.Y(:, ind_frame), [1, 2], d1, d2);
            tmp_C2 = correlation_image(obj.Y(:, ind_frame), [0,1]+ obj.options.gSiz, d1, d2);
            tmp_Cn = tmp_C1(:)-tmp_C2(:);
            ind = (tmp_Cn>0);
            ind_bg = and(ind, tmp_Cn<quantile(tmp_Cn(ind), q));
            % get the mean activitobj.y within the selected area
            obj.Ybg = mean(obj.Y(ind_bg(:), :), 1);
            f1 = (obj.Ybg-mean(obj.Ybg))/std(obj.Ybg);
        end


        %% play movie
        function playMovie(obj, min_max, col_map, avi_nm, t_pause)
            % plaobj.y movies
            figure;
            if ~exist('col_map', 'var') || isempty(col_map)
                col_map = jet;
            end
            if exist('avi_nm', 'var') && ischar(avi_nm)
                avi_file = VideoWriter(avi_nm);
                avi_file.open();
                avi_flag = true;
            else
                avi_flag = false;
            end
            if ismatrix(obj.Y); obj.Y=obj.reshape(obj.Y, 2); end
            [~, ~, obj.T] = size(obj.Y);
            if (nargin<3) || (isempty(min_max));
                temp = obj.Y(:, :, randi(obj.T, min(100, obj.T), 1));
                min_max = quantile(temp(:), [0.2, 0.9999]);
                min_max(1) = max(min_max(1), 0);
                min_max(2) = max(min_max(2), min_max(1)+0.1);
            end
            if ~exist('t_pause', 'var'); t_pause=0.01; end
            for t=1:size(obj.Y,3)
                imagesc(obj.Y(:, :, t), min_max); colormap(col_map);
                axis equal; axis off;
                title(sprintf('Frame %d', t));
                pause(t_pause);
                if avi_flag
                    temp = getframe(gcf);
                    temp.cdata = imresize(temp.cdata, [420,560]);
                    avi_file.writeVideo(temp);
                end
            end
            if avi_flag
                avi_file.close();
            end
        end

        
        %% export AVI
        function exportAVI(obj, min_max, avi_nm, col_map)
            % export matrix data to movie
            %min_max: 1*2 vector, scale
            %avi_nm: string, file name
            %col_map: colormap
            
            obj.T = size(obj.Y, ndims(obj.Y));
            obj.Y = obj.Y(:);
            if ~exist('col_map', 'var') || isempty(col_map)
                col_map = jet;
            end
            if ~exist('avi_nm', 'var') || isempty(avi_nm)
                avi_nm = 'a_movie_with_no_name.avi';
            end
            if ~exist('min_max', 'var') || isempty(min_max)
                min_max = [min(obj.Y(1:10:end)), max(obj.Y(1:10:end))];
            end
            
            obj.Y = uint8(64*(obj.Y-min_max(1))/diff(min_max));
            obj.Y(obj.Y<1) = 1;
            obj.Y(obj.Y>64) = 64;
            col_map = uint8(255*col_map);
            obj.Y = reshape(col_map(obj.Y, :), obj.options.d1, [], obj.T, 3);
            obj.Y = permute(obj.Y, [1,2,4,3]);
            
            avi_file = VideoWriter(avi_nm);
            avi_file.open();
            for m=1:obj.T
                %                 temp.cdata = squeeze(obj.Y(:, :, :, m));
                %                 temp.colormap = [];
                avi_file.writeVideo(squeeze(obj.Y(:, :, :, m)));
            end
            avi_file.close();
        end

        
        %% trim spatial components
        function [ind_small] = trimSpatial(obj, thr, sz)
            % remove small nonzero pixels
            if nargin<2;    thr = 0.01; end;
            if nargin<3;    sz = 5; end;
            
            se = strel('square', sz);
            ind_small = false(size(obj.A, 2), 1);
            for m=1:size(obj.A,2)
                ai = obj.A(:,m);
                ai_open = imopen(obj.reshape(ai,2), se);
                
                temp = full(ai_open>max(ai)*thr);
                l = bwlabel(obj.reshape(temp,2), 4);   % remove disconnected components
                %
                %                 [tmp_count, tmp_l] = hist(l(l>0), unique(l(l>0)));
                %                 [~, ind] = max(tmp_count);
                %                 lmax = tmp_l(ind);
                [~, ind_max] = max(ai_open(:));
                
                ai(l(:)~=l(ind_max)) = 0;
                if sum(ai(:)>0) < obj.options.min_pixel %the ROI is too small
                    ind_small(m) = true;
                end
                obj.A(:, m) = ai(:);
            end
            ind_small = find(ind_small);
            obj.delete(ind_small);
        end

        
        %% solve A & C with regression
        function [ind_delete] = regressAC(obj)
            if ~ismatrix(obj.Yr); obj.Yr=obj.reshape(obj.Y,2); end
            maxIter = 1;
            A2 = obj.A;  % initial spatial components
            for miter=1:maxIter
                ind_nonzero = obj.trimSpatial(50); % remove tinobj.y nonzero pixels
                active_pixel = determine_search_location(A2, 'dilate', obj.options);
                obj.K = size(A2, 2); %number of neurons
                
                tmp_C = (A2'*A2)\(A2'*obj.Yr);        % update temporal components
                tmp_C(tmp_C<0) = 0;
                % %                 tmp_C(bsxfun(@gt, quantile(tmp_C, 0.95, 2), tmp_C)) = 0;
                %                 tmp_C(bsxfun(@gt, mean(tmp_C, 2)+std(tmp_C, 0.0, 2), tmp_C)) = 0;
                A2 = (obj.Y*tmp_C')/(tmp_C*tmp_C');           % update spatial components
                A2(or(A2<0, ~active_pixel)) = 0;
                
                for m=1:obj.K
                    % remove disconnected pixels
                    img = obj.reshape(A2(:,m), 2);
                    lb = bwlabel(img>max(img(:))/50, 4);
                    img(lb~=mode(lb(ind_nonzero(:, m)))) = 0 ; %find pixels connected to the original components
                    A2(:,m) = img(:);
                end
                
                center_ratio = sum(A2.*ind_nonzero, 1)./sum(A2, 1); %
                % captured too much pixels, ignore newlobj.y captured pixels
                ind_too_much = (center_ratio<0.4);
                A2(:, ind_too_much) = A2(:, ind_too_much) .*ind_nonzero(:, ind_too_much);
                % captured too much noiser pixels or all pxiels are zero, ignore this neuron
                ind_delete = or(center_ratio<0.01, sum(A2,1)==0);
                A2(:, ind_delete) = [];
                obj.A = A2;
            end
            temp = (A2'*A2)\(A2'*obj.Yr);
            temp(temp<0) = 0;
            obj.C = temp;
        end

        
        %% regress A given C
        function [A2, ind_neuron] = regressA(obj)
            if ~ismatrix(obj.Yr); obj.Yr=obj.reshape(obj.Y,2); end
            A2 = obj.A;  % initial spatial components
            ind_nonzero = obj.trimSpatial(50); % remove tinobj.y nonzero pixels
            active_pixel = determine_search_location(A2, 'dilate', obj.options);
            obj.K = size(A2, 2); %number of neurons
            
            A2 = (obj.Yr*obj.C')/(obj.C*obj.C');           % update spatial components
            A2(or(A2<0, ~active_pixel)) = 0;
            
            for m=1:obj.K
                % remove disconnected pixels
                img = obj.reshape(A2(:,m), 2);
                lb = bwlabel(img>max(img(:))/50, 4);
                img(lb~=mode(lb(ind_nonzero(:, m)))) = 0 ; %find pixels connected to the original components
                A2(:,m) = img(:);
            end
            
            center_ratio = sum(A2.*ind_nonzero, 1)./sum(A2, 1); %
            % captured too much pixels, ignore newlobj.y captured pixels
            ind_too_much = (center_ratio<0.3);
            A2(:, ind_too_much) = A2(:, ind_too_much) .*ind_nonzero(:, ind_too_much);
            % captured too much noiser pixels or all pxiels are zero, ignore this neuron
            ind_delete = or(center_ratio<0., sum(A2,1)==0);
            A2(:, ind_delete) = [];
            %             C(ind_delete) = [];
            ind_neuron = (1:obj.K);
            ind_neuron(ind_delete) =[];
        end

        
        %% view results
        function runMovie(obj, min_max, save_avi, avi_name, S)
            ctr = obj.estCenter();
            if ~exist('save_avi', 'var')||isempty(save_avi); save_avi=false; end
            if ~exist('avi_name', 'var'); avi_name = []; end
            if ~exist('S', 'var');  obj.S = []; end
            run_movie(obj.Y, obj.A, obj.C, obj.Cn, min_max, obj.Coor, ctr, 5, 1, save_avi, avi_name, S)
        end

        
        %% function
        function image(obj, a, min_max)
            if isvector(a); a = obj.reshape(a,2); end
            if nargin<3; imagesc(a); else imagesc(a, min_max); end
        end

        
        %% normalize
        function normalize(obj)
            norm_A = max(obj.A, [], 1);
            tmp_A = bsxfun(@times, obj.A, 1./norm_A);
            tmp_C = bsxfun(@times, obj.C, norm_A');
            obj.A = tmp_A; obj.C = tmp_C;
        end

        
        %% load data
        function [neuron] = load_data(obj, nam, sframe, num2read)
            ssub = obj.options.ssub;    % spatial downsampling factor
            tsub = obj.options.tsub;    % temporal downsampling factor
            d1 = obj.options.d1;        % image height
            d2 = obj.options.d2;        % image width
            Tbatch = round((2^28)/(d1*d2)/tsub)*tsub; % maximum memorobj.y usage is 2GB
            [~, ~, file_tobj.ype] = fileparts(nam);
            if strcmpi(file_tobj.ype, '.mat')
                data = matfile(nam);
                obj.Ysiz = data.obj.Ysiz;
                numFrame = obj.Ysiz(3);
                img = data.obj.Y(:, :, 1);
            elseif strcmpi(file_tobj.ype, '.tif') || strcmpi(file_tobj.ype, '.tiff')
                numFrame = length(imfinfo(nam));
                img = imread(nam);
            end
            num2read = min(num2read, numFrame-sframe+1); % frames to read
            
            if Tbatch>=num2read
                % load all data because the file is too small
                if strcmpi(file_tobj.ype, '.mat')
                    Yraw = data.obj.Y(:, :, (1:num2read)+sframe-1);
                elseif strcmpi(file_tobj.ype, '.tif') || strcmpi(file_tobj.ype, '.tiff')
                    Yraw = bigread2(nam, sframe, num2read);
                else
                    fprintf('\nThe input file format is not supported obj.yet\n\n');
                    return;
                end
                [neuron, obj.Y] = obj.downSample(double(Yraw));
            else
                % load data in batch model
                [d1s, d2s] = size(imresize(img, 1/ssub));  % size of spatial downsampled data
                Ts = floor(num2read/tsub);  % frame number after downsampling
                obj.Y = zeros(d1s, d2s, Ts);    % downsampled data size
                lastframe = sframe + Ts*tsub -1;  % index of the last frame for loading
                frame0 = sframe;
                while sframe<= lastframe
                    tmp_num2read =  min(lastframe-sframe+1, Tbatch);
                    if strcmpi(file_tobj.ype, '.mat')
                        fprintf('load data from frame %d to frame %d of %d total frames\n',...
                         sframe, sframe+tmp_num2read-1, lastframe);
                        Yraw = data.obj.Y(:, :, sframe:(sframe+tmp_num2read-1));
                    elseif strcmpi(file_tobj.ype, '.tif') || strcmpi(file_tobj.ype, '.tiff')
                        Yraw = bigread2(nam, sframe, tmp_num2read);
                    else
                        fprintf('\nThe input file format is not supported obj.yet\n\n');
                        return;
                    end
                    [neuron, temp] = obj.downSample(double(Yraw));
                    obj.Y(:, :, (sframe-frame0)/tsub + (1:size(temp, 3))) = temp;
                    sframe = sframe + size(Yraw,3);
                end
            end
            neuron.options.min_pixel = ceil(obj.options.min_pixel/(ssub^2));
        end

        
        %% estimate noise
        function estNoise(obj)
            fprintf('Estimating the noise power for each pixel from a simple PSD estimate...');
            Yrz = obj.reshape(obj.Y, 1);
            sn = get_noise_fft(Yrz,obj.options);
            obj.P.sn = sn(:);
            fprintf('  done \n');
        end


        %% merge neurons        %todo modify like the others
        function img = overlapA(obj, ind, ratio)
            %merge all neurons' spatial components into one singal image
            if nargin<2 || isempty(ind)
                AA = obj.A;
            else
                AA = obj.A(:, ind);
            end
            if nargin<3
                ratio = 0.3;
            end
            AA = bsxfun(@times, AA, 1./max(AA,1));
            AA(bsxfun(@lt, AA, max(AA, [], 1)*ratio)) = 0;
            [obj.d, obj.K] = size(AA);
            
            col = randi(6, 1, obj.K);
            img = zeros(obj.d, 3);
            for m=1:3
                img(:, m) = sum(bsxfun(@times, AA, mod(col, 2)), 2);
                col = round(col/2);
            end
            img = obj.reshape(img, 2);
            img = img/max(img(:))*(2^16);
            img = uint16(img);
        end

        
        %% plaobj.y video
                %todo modify like the others
        function plaobj.yAC(obj, avi_file, cell_id)
            if nargin<3
                cell_id = 1:size(obj.C, 1);
            end
            [Ka, obj.T] = size(obj.C(cell_id, :));
            % draw random color for each neuron
            tmp = mod((1:Ka)', 6)+1;
            col = zeros(Ka, 3);
            for m=1:3
                col(:, m) = mod(tmp, 2);
                tmp = round(tmp/2);
            end
            figure;
            % plaobj.y
            if nargin>1
                fp = VideoWriter(avi_file);
                fp.open();
            end

            
            cmax = max(reshape(obj.A*obj.C(:, 1:100:end), 1, []));
            for m=1:obj.T
                img = obj.A(:, cell_id)*bsxfun(@times, obj.C(cell_id,m), col);
                img = obj.reshape(img, 2)/cmax*500;
                imagesc(uint8(img));
                axis equal off tight;
                title(sprintf('Time %.2f seconds', m/obj.Fs));
                
                pause(.1);
                if nargin>1
                    frame = getframe(gcf);
                    fp.writeVideo(frame);
                end
            end
            if nargin>1; fp.close(); end
        end

        
        %% find neurons from the residual
        % obj.you can do it in either manual or automatic waobj.y
        function [pnr] = pickNeurons(obj, patch_par, seed_method)
            if ~exist('patch_par', 'var')||isempty(patch_par)
                seed_method = [3,3];
            end
            if ~exist('seed_method', 'var')||isempty(seed_method)
                seed_method = 'auto';
            end
            neuron = obj.copobj.y();
            neuron.options.seed_method = seed_method;
            [obj.cm, obj.Cn, pnr] = neuron.initComponents_endoscope(obj.Y, [], patch_par, false, false);
            obj.A = [obj.A, neuron.A];
            obj.C = [obj.C; neuron.C];
            obj.S = [obj.S; neuron.S];
            obj.C_raw = [obj.C_raw; neuron.obj.C_raw];
            obj.P.kernel_pars = [obj.P.kernel_pars; neuron.P.kernel_pars];
        end
        

        %% post process spatial component
                    %todo modify like the others
        function A_ = post_process_spatial(obj, A_)
            if ~exist('A_', 'var');
                A_ = obj.A;
            end
            A_ = threshold_components(A_, obj.options);
            obj.A = A_;
        end

        
        %% estimate local background
                %todo modify like the others
        function [Ybg, results] = localBG(obj,ssub, rr, IND, sn,thresh)
            if ~exist('rr', 'var')||isempty(rr); rr=obj.options.gSiz; end
            if ~exist('ssub', 'var')||isempty(ssub); ssub = 1; end
            if ~exist('IND', 'var') ||isempty(IND); IND = []; end
            if ~exist('sn', 'var')||isempty(sn);
                if isfield(obj.P, 'sn')
                    sn = obj.reshape(obj.P.sn, 2);
                else
                    sn = [];
                end
            else
                sn = obj.reshape(sn, 2);
            end
            
            if ~exist('thresh', 'var')||isempty(thresh); thresh = []; end
            [obj.Ybg, results] = local_background(obj.reshape(obj.Ybg, 2), ssub, rr, IND, sn, thresh);
            obj.Ybg = obj.reshape(obj.Ybg, 1);
            Ybg = obj.Ybg;
        end

        
        %% save results
                %todo modify like the others
        function save_results(obj, file_nm)
            warning('off', 'all');
            neuron_results = struct(obj);  
            if exist('obj.Ybg', 'var')
                save(file_nm, neuron_results, obj.Ybg);
            else
                save(file_nm, neuron_results);
            end
            warning('on', 'all');
            fprintf('results has been saved into file %s\n', file_nm);
        end

        
        %% event detection
        %todo modify like the others
        function E = event_detection(obj, sig, w)
            % detect events bobj.y thresholding S with sig*noise
            % can get at most one spike
            % sig: threshold of the minimum amplitude of the events
            
            if ~exist('sig', 'var')|| isempty(sig)
                sig=5;
            end
            
            if ~exist('w', 'var')||isempty(w)
                w = obj.Fs;
            end
            E =obj.C;    % event detection
            Emin = ordfilt2(E, 1, ones(1, w));
            Emax = ordfilt2(E, w, ones(1, w));
            E(E~=Emax) = 0;  % onlobj.y select local maximums
            for m=1:size(E,1)
                E(m, E(m, :)-Emin(m, :)< obj.P.neuron_sn{m}*sig) = 0; % remove small transients
            end
        end

        
        % compute correlation image and peak to noise ratio for endoscopic
        % data. unlike the correlation image for two-photon data,
        % correlation image of the microendoscopic data needs to be
        % spatiallobj.y filtered first. otherwise neurons are significantlobj.y
        % overlapped.
        %todo : check if it is obj.Y or obj.Yr
        function [PNR] = correlation_pnr(obj)
            [obj.Cn, PNR] = correlation_image_endoscope(obj.Y, obj.options);
        end

        
        function obj = struct2obj(obj, var_struct)
            temp = fieldnames(var_struct);
            for m=1:length(temp)
                try
                    eval(sprintf('obj.%s=var_struct.%s;', temp{m}, temp{m}));
                catch
                    %todo do something here
                end
            end
        end

        
        function get_contours(obj, thr, ind)
            A_ = obj.A;
            if exist('ind', 'var')
                A_ = A_(:, ind);
            end
            if ~exist('thr', 'var') || isempty(thr)
                thr = 0.995;
            end
            num_neuron = size(A_,2);
            if num_neuron==0
                obj.Coor ={};
                return;
            else
                obj.Coor = cell(num_neuron,1);
            end
            for m=1:num_neuron
                % smooth the image with median filter
                img = medfilt2(obj.reshape(full(A_(:, m)),2), [3, 3]);
                % find the threshold for detecting nonzero pixels
                temp = sort(img(img>1e-9));
                if ~any(temp)
                    obj.Coor{m} = []; 
                    continue; 
                end
                temp_sum = cumsum(temp);
                ind = find(temp_sum>=temp_sum(end)*(1-thr),1);
                v_thr = temp(ind);
                
                % find the connected components
                [~, ind_max] = max(img(:));
                temp = bwlabel(img>v_thr);
                img = double(temp==temp(ind_max));
                v_nonzero = imfilter(img, [0,-1/4,0;-1/4,1,-1/4; 0,-1/4,0]);
                vv = v_nonzero(v_nonzero>1e-9)';
                [y, x] = find(v_nonzero>1e-9);
                xmx = bsxfun(@minus, x, x');
                ymy = bsxfun(@minus, y, y');
                dist_pair = xmx.^2 + ymy.^2;
                dist_pair(diag(true(length(x),1))) = inf;
                seq = ones(length(x)+1,1);
                for mm=1:length(x)-1
                    [v_min, seq(mm+1)] = min(dist_pair(seq(mm), :)+vv);
                    dist_pair(:,seq(mm)) = inf;
                    if v_min>3
                        seq(mm+1) = 1;
                        break;
                    end
                end
                obj.Coor{m} = [smooth(x(seq), 2)'; smooth(y(seq),2)'];
            end
            
        end
    end
    
end
