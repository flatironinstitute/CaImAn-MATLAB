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
    % UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
    'search_method      ' % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'dilate')
    'spatial_parallel   ' % update pixels in parallel (default: 1 if present)
    % determine_search_location.m
    'min_size           ' % minimum size of ellipse axis (default: 3)
    'max_size           ' % maximum size of ellipse axis (default: 8)
    'dist               ' % expansion factor of ellipse (default: 3)
    'se                 ' % morphological element for dilation (default: strel('disk',4,0))
    % threshold_components.m
    'nrgthr             ' % energy threshold (default: 0.9999)
    'clos_op            ' % morphological element for closing (default: strel('square',3))
    'medw               ' % size of median filter (default: [3,3])
    % UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
    'deconv_method      '    % method for spike deconvolution (default: 'constrained_foopsi')
    'restimate_g        '    % flag for updating the time constants for each component (default: 1)
    'temporal_iter      '    % number of block-coordinate descent iterations (default: 2)
    'temporal_parallel  ' % flag for parallel updating of temporal components (default: true if present)
    'full_A             ' % if true turn A into full matrix. If false turn Y into double precision (default: false)
    % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
    'method             ' % methods for performing spike inference ('dual','cvx','spgl1','lars') (default:'cvx')
    'bas_nonneg         ' % flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y) (default 1)
    'fudge_factor       ' % scaling constant to reduce bias in the time constant estimation (default 1 - no scaling)
    'resparse           ' % number of times that the solution is resparsened (default: 0)
    % MERGING (merge_ROIs.m)
    'merge_thr          ' % merging threshold (default: 0.85)
    'fast_merge         ' % flag for using fast merging (default 1)
    % DF/F (extract_DF_F.m)
    'df_prctile         ' % percentile to be defined as baseline (default 50, median)
    'df_window          ' % length of running window (default [], no window)
    % CONTOUR PLOTS (plot_contours.m)
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
    % CLASSIFY COMPONENTS (classify components.m)
    'cl_thr             ' % overlap threshold for energy for a component to be classified as true (default: 0.8)
    % ORDER COMPONENTS (order_components.m)
    'nsd                ' % number of standard deviations (default: 3)
    'nfr                ' % number of consecutive frames (default: 3)
    ];

[m,n] = size(Names);
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
    {'logmexp'}
    {3000}
    {false}
    {5}
    {false}
    {[]}
    {false}
    {[64,64]}
    {false}
    % UPDATING SPATIAL COMPONENTS (unpdate_spatial_components.m)
    {'dilate'}
    {~isempty(which('parpool'))}
    % determine_search_location.m
    {3}
    {8}
    {3}
    {strel('disk',4,0)}
    % threshold_components.m
    {0.995}
    {strel('square',3)}
    {[3,3]}
    % UPDATING TEMPORAL COMPONENTS (update_temporal_components.m)
    {'constrained_foopsi'}
    {1}
    {2}
    {~isempty(which('parpool'))}
    {false}
    % CONSTRAINED DECONVOLUTION (constrained_foopsi.m)
    {'cvx'}
    {1}
    {0.99}
    {0}
    % MERGING (merge_ROIs.m)
    {0.85}
    {1}
    % DF/F (extract_DF_F.m)
    {50}
    {[]}
    % CONTOUR PLOTS (plot_contours.m)
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
    % CLASSIFY COMPONENTS (classify_components.m)
    {0.8}
    % ORDER COMPONENTS (order_components.m)
    {3}
    {5}
    ];

for j = 1:m
    if eval(['isempty(options.' Names(j,:) ')'])
        eval(['options.' Names(j,:) '= Values{j};']);
    end
end
