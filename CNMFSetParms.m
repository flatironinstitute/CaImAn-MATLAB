function options = CNMFSetParms(varargin)

% Struct for setting the CNMF algorithm parameters. Any parameter that is
% not set gets a default value

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2015

Names = [
    % dataset info
    'd1                 ' % number of rows
    'd2                 ' % number of cols
    'd3                 ' % number of planes (for 3d imaging, default: 1)
    % INITIALIZATION  (initialize_components.m)
    'ssub               ' % spatial downsampling factor (default: 1)
    'tsub               ' % temporal downsampling factor (default: 1)
    'init_method        ' % initialization method ('greedy','greedy_corr','sparse_NMF','HALS') (default: 'greedy')
    'rem_prct           ' % percentile to be removed before initialization (default: 20)
    'noise_norm         ' % normalization by noise estimate prior to initialization (default: true)
    'noise_norm_prctile ' % minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
    % greedy_corr parameters (greedyROI_corr.m)
    'min_corr           ' % minimum local correlation for initializing a neuron (default: 0.3)
    % greedyROI parameters (greedyROI.m)
    'gSig               ' % half size of neurons to be found (default: [5,5])
    'gSiz               ' % half size of bounding box for each neuron (default: 2*gSig+1)
    'nb                 ' % number of background components (default: 1)
    'nIter              ' % maximum number of rank-1 NMF iterations during refining
    'med_app            ' % number of timesteps to be interleaved for fast (approximate) median calculation (default: 1)
    'save_memory        ' % process data sequentially to save memory (default: 0)
    'chunkSiz           ' % filter this number of timesteps each time (default: 100)
    'windowSiz          ' % size of window over which is computed sequentially (default: 32 x 32)
    'rolling_sum        ' % flag for using rolling sum to detect new components (default: True)
    'rolling_length     ' % length of rolling window (default: 100)
    % sparse_NMF parameters (sparse_NMF_initialization.m)
    'snmf_max_iter      ' % max # of sparse NMF iterations
    'err_thr            ' % relative change threshold for stopping sparse_NMF
    'eta                ' % frobenious norm factor *max(Y(:))^2
    'beta               ' % sparsity factor
    % HALS initialization parameters (HALS_initialization.m)
    'max_iter_hals_in   ' % maximum number of HALS iterations
    % HALS parameters (HALS_2d.m)
    'bSiz               ' % expand kernel for HALS growing (default: 3)
    'maxIter            ' % maximum number of HALS iterations (default: 5)
    % Noise and AR coefficients calculation (preprocess_data.m)
    'noise_range        ' % frequency range over which to estimate the noise (default: [0.25,0.5])
    'noise_method       ' % method for which to estimate the noise level (default: 'logmexp')
    'max_timesteps      ' % maximum number of timesteps over which to estimate noise (default: 3000)
    'flag_g             ' % compute global AR coefficients (default: false)
    'lags               ' % number of extra lags when computing the AR coefficients (default: 5)
    'include_noise      ' % include early lags when computing AR coefs (default: 0)
    'pixels             ' % pixels to include when computing the AR coefs (default: 1:numel(Y)/size(Y,ndims(Y)))
    'split_data         ' % split data into patches for memory reasons (default: 0)
    'block_size         ' % block size for estimating noise std in patches (default: [64,64])
    'cluster_pixels     ' % cluster pixels to active/inactive based on the PSD density (default: false)
    'extract_max        ' % extract the maximum activity intervals for each pixel (default: false)
    'max_nlocs          ' % number of local maxima to be extracted (default: 10)
    'max_width          ' % length of each interval (default: 11)
    % UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
    'spatial_method     ' % method for updating spatial components 'constrained' or 'regularized' (default: 'regularized')
    'search_method      ' % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'dilate')
    'spatial_parallel   ' % update pixels in parallel (default: 1 if present)
    % determine_search_location.m
    'min_size           ' % minimum size of ellipse axis (default: 3)
    'max_size           ' % maximum size of ellipse axis (default: 8)
    'dist               ' % expansion factor of ellipse (default: 3)
    'se                 ' % morphological element for dilation (default: strel('disk',1,0))
    % threshold_components.m
    'thr_method         ' % method to threshold ('max' or 'nrg', default 'max')
    'maxthr             ' % threshold of max value below which values are discarded (default: 0.25)
    'nrgthr             ' % energy threshold (default: 0.995)
    'clos_op            ' % morphological element for closing (default: strel('square',3))
    'medw               ' % size of median filter (default: [3,3])
    'conn_comp          ' % extract largest connected component (binary, default: true)
    % UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
    'p                  ' % order of AR model dynamics (default: 1)
    'deconv_method      ' % method for spike deconvolution (default: 'constrained_foopsi')
    'restimate_g        ' % flag for updating the time constants for each component (default: 1)
    'temporal_iter      ' % number of block-coordinate descent iterations (default: 2)
    'temporal_parallel  ' % flag for parallel updating of temporal components (default: true if present)
    'full_A             ' % if true turn A into full matrix. If false turn Y into double precision (default: false)
    % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
    'method             ' % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
    'bas_nonneg         ' % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
    'fudge_factor       ' % scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)
    'resparse           ' % number of times that the solution is resparsened (default: 0)
    % OASIS penalty parameters (deconvolveCa.m)
    'lam_pr             ' % false positive probability for determing lambda penalty
    'spk_SNR            ' % spike SNR for min spike value
    % MERGING (merge_ROIs.m)
    'merge_thr          ' % merging threshold (default: 0.85)
    'fast_merge         ' % flag for using fast merging (default 1)
    % DF/F (extract_DF_F.m)
    'df_prctile         ' % percentile to be defined as baseline (default 20)
    'df_window          ' % length of running window (default [], no window)
    % CONTOUR PLOTS (plot_contours.m)
    'plot_bck_image     ' % plot background image or overlay on existing one (deafult: true)
    'cont_threshold     ' 
    % VIDEO (make_patch_video.m)
    'ind                ' % indeces of components to be shown (deafult: 1:4)
    'skip_frame         ' % skip frames when showing the video (default: 1 (no skipping))
    'sx                 ' % half size of representative patches (default: 16)
    'make_avi           ' % flag for saving avi video (default: 0)
    'show_background    ' % flag for displaying the background in the denoised panel (default: 1)
    'show_contours      ' % flag for showing the contour plots of the patches in the FoV (default: 0)
    'cmap               ' % colormap for plotting (default: 'default')
    'name               ' % name of saved video file (default: based on current date)
    % PLOT COMPONENTS (view_patches.m)
    'plot_df            ' % flag for displaying DF/F estimates (default: 1)
    'make_gif           ' % save animation (default: 0)
    'save_avi           ' % save video (default: 0)
    'pause_time         ' % time to pause between each component (default: Inf, user has to click)
    % CLASSIFY COMPONENTS PIXELS (classify_components_pixels.m)
    'cl_thr             ' % overlap threshold for energy for a component to be classified as true (default: 0.8)
    % CLASSIFY COMPONENTS with CORRELATION (classify_comp_corr.m)
    'space_thresh       ' % threshold for r-value in space (default: 0.4)
    'time_thresh        ' % threshold for r-value in time (default: 0.4)
    'A_thresh           ' % threshold for determining overlap (default: 0.1)
    'Npeaks             ' % # of peaks to be considered (default: 5)
    'peak_int           ' % interval around the peak (default: -2:5)
    'MinPeakDist        ' % minimum peak distance for finding points of high activity  (default: 10)
    % ORDER COMPONENTS (order_components.m)
    'nsd                ' % number of standard deviations (default: 3)
    'nfr                ' % number of consecutive frames (default: 3)
    % PATCHES          (run_CNMF_patches.m)
    'gnb                ' % number of global background components (default: 1)
    'create_memmap      ' % create a memory mapped file if it is not provided in the input (default: false)    
    'classify_comp      ' % classify components based on correlation values (default: true)
    'refine_flag        ' % refine components within patch processing after merging (default: true)    
    'patch_space_thresh ' % space correlation threshold within patch (default: 0.2)
    'patch_time_thresh  ' % time correlation threshold within patch (default: 0.25)
    'patch_min_SNR      ' % minimum SNR for accepting exceptional events within patch (default: 0.5)
    'patch_min_fitness  ' % maximum fitness threshold within patch (default: log(normcdf(-patch_min_SNR))*N_samples_exc)
    'patch_min_fit_delta' % maximum fitness_delta threshold within patch (default: -2)
    'patch_cnn_thr      ' % threshold for CNN classifier within a patch (default: 0.05)
    % parameters for microendoscope 
    'min_pnr            '
    'seed_method        '    
    'min_pixel          ' % minimum number of nonzero pixels for a neuron 
    'bd                 ' % number of pixels to be ignored in the boundary 
    'deconv_flag        ' % perform deconvolution or not     
    % parameters for max probability test (trace_fit_extreme.m)
    'max_pr_thr         ' % threshold for keeping components (default: 0.9)
    'fr                 ' % imaging frame rate in Hz (defaut: 30)
    't_int              ' % length of each trial in sec (default: 0.25)
    'sn_fac             ' % multiplicative factor for estimated noise level (default: 1)
    % parameters for thresholding based on size (classify_components.m)
    'max_size_thr       ' % maximum size of each component in pixels (default: 300)
    'min_size_thr       ' % minimum size of each component in pixels (default: 9)
    'size_thr           ' % fraction of max value for thresholding each component before determining its size (default 0.2)
    % parameters for registering components across different sessions (register_ROIs.m)
    'dist_exp           ' % exponent for calculating the distance between different ROIs (default: 1)
    'dist_thr           ' % distance threshold above which dist = Inf (default: 0.5)
    'dist_maxthr        ' % max thresholding for components before turing into binary masks (default: 0.15)
    'dist_overlap_thr   ' % threshold for detecting if one ROI is a subset of another (deafult: 0.8)
    'plot_reg           ' % plot registered ROIs (default: true)
    % parameters for computing event exceptionality (compute_event_exceptionality.m)
    'min_SNR            ' % minimum SNR for accepting exceptional events (default: 2)
    'decay_time         ' % length of a typical transient in seconds
    'robust_std         ' % use robust std for computing noise in traces (false)
    'N_samples_exc      ' % number of samples over which to compute (default: ceil(decay_time*fr))
    'min_fitness        ' % threshold on time variability  (default: log(normcdf(-min_SNR))*N_samples_exc)    
    'min_fitness_delta  ' % threshold on the derivative of time variability
    % parameters for CNN classifier (cnn_classifier.m)
    'cnn_thr            ' % threshold for CNN classifier (default: 0.2)
    ];

[m,~] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if ischar(arg), break; end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(sprintf(['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with OPTIMSET.'], i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
                eval(['val = arg.' Names(j,:) ';']);
            else
                val = [];
            end
            if ~isempty(val)
                eval(['options.' Names(j,:) '= val;']);
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string parameter name.', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized parameter name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' deblank(Names(j(1),:))];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names(k,:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;                      % we expect a value next
        
    else
        eval(['options.' Names(j,:) '= arg;']);
        expectval = 0;
        
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end

Values = [
    % dataset info
    {[]}
    {[]}
    {1}
    % INITIALIZATION  (initialize_components.m)
    {1}
    {1}
    {'greedy'}
    {20}
    {true}
    {2}
    % greedy_corr parameters (greedyROI_corr.m)
    {.3}
    % greedyROI parameters (greedyROI.m)
    {5}
    {[]}
    {1}
    {5}
    {1}
    {0}
    {100}
    {[32,32]}
    {true}
    {100}
    % sparse_NMF parameters (sparse_NMF_initialization.m)
    {100}
    {1e-4}
    {1}
    {.5}
    % HALS initialization parameters (HALS_initialization.m)
    {5}
    % HALS parameters (HALS_2d.m)
    {3}
    {5}
    % Noise and AR coefficients calculation (preprocess_data.m)
    {[0.25,0.5]}
    {'mean'}
    {3000}
    {false}
    {5}
    {false}
    {[]}
    {false}
    {[64,64]}
    {false}
    {false}
    {30}
    {21}
    % UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
    {'regularized'}
    {'dilate'}
    {~isempty(which('parpool'))}
    % determine_search_location.m
    {3}
    {8}
    {3}
    {strel('disk',1,0)}
    % threshold_components.m
    {'max'}
    {0.25}
    {0.995}
    {strel('square',3)}
    {[3,3]}
    {true}
    % UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
    {1}
    {'constrained_foopsi'}
    {1}
    {4}
    {~isempty(which('parpool'))}
    {false}
    % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
    {'cvx'}
    {1}
    {0.98}
    {0}
    % OASIS penalty parameters (deconvolveCa.m)
    {0.99}
    {0.5}
    % MERGING (merge_ROIs.m)
    {0.85}
    {1}
    % DF/F (extract_DF_F.m)
    {20}
    {[]}
    % CONTOUR PLOTS (plot_contours.m)
    {true}
    {0.9}
    % VIDEO (make_patch_video.m)
    {1:4}
    {1}
    {16}
    {0}
    {1}
    {1}
    {'default'}
    {['video_',datestr(now,30),'.avi']}
    % PLOT COMPONENTS (plot_patches.m)
    {1}
    {0}
    {0}
    {Inf}
    % CLASSIFY COMPONENTS PIXELS (classify_components_pixels.m)
    {0.8}
    % CLASSIFY COMPONENTS with CORRELATION (classify_comp_corr.m)
    {0.4}
    {0.4}
    {0.1}
    {5}
    {-2:6}
    {10}
    % ORDER COMPONENTS (order_components.m)
    {3}
    {5}
    % PATCHES          (run_CNMF_patches.m)
    {1}
    {false}    
    {true}
    {true}
    {0.2}
    {0.25}
    {0.5}
    {[]}
    {-2}
    {0.05}
    % parameters for microendoscope
    {10}
    {'auto'}
    {5}
    {3}    
    {true}
    % parameters for max probability test (trace_fit_extreme.m)
    {0.9}
    {30}
    {0.25}
    {1}
    % parameters for size based thresholding (classify_components.m)
    {320}
    {9}
    {0.2}
    % parameters for registering components across different sessions (register_ROIs.m)
    {1}
    {0.5}
    {0.15}
    {0.8}
    {true}
    % parameters for computing event exceptionality (compute_event_exceptionality.m)
    {2}
    {0.4}
    {false}
    {[]}
    {[]}
    {-5}
    % parameters for CNN classifier (cnn_classifier.m)
    {0.2}
    ];

for j = 1:m
    if eval(['isempty(options.' Names(j,:) ')'])
        eval(['options.' Names(j,:) '= Values{j};']);
    end
end

if isempty(options.N_samples_exc); options.N_samples_exc = ceil(options.fr*options.decay_time); end
if isempty(options.min_fitness); options.min_fitness = log(normcdf(-options.min_SNR))*options.N_samples_exc; end
if isempty(options.patch_min_fitness); options.patch_min_fitness = log(normcdf(-options.patch_min_SNR))*options.N_samples_exc; end
