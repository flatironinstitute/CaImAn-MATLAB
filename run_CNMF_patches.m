function [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options)
% RUN_CNMF_PATCHES - apply CNMF algorithm on overlapping patches in parallel
%
%   [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options)
%
% Run the constrained NMF algorithm on a large dataset by operating on
% spatially overlapping patches in parallel and then merging the results.
% The inputs is memory mapped, allowing for large datasets to be processed
% with reduced memory requirements. Processing in patches also allows the
% identification of weaker neurons without the need of normalization.
% The components are also classified by retaining only the components that
% correlate well with the raw data through classify_comp_corr.m
%
% INPUTS:
% data:    .mat file containing
%            data.Y      (the data matrix in the original dimensions)
%            data.Yr     (the data matrix reshaped in 2d format)
%            data.sizY   (dimensions of the original dataset)
%            data.nY     (minimum value of dataset)
%          OR the original dataset in 3d/4d format in which case the user
%          chooses whether to create a memory mapped file
% K:       number of components to be found in each patch
% patches: cell array containing the start and end points of each patch
% tau:     half-size of each cell for initializing the components
% p:       order of autoregressive progress
% options: struct for algorithm parameters
%
% OUTPUTS:
% A:       Matrix of spatial components
% b:       Spatial background
% C:       Matrix of temporal components
% f:       Temporal background
% P:       Struct for neuron parameters
% RESULTS: Results of the CNMF algorithm on individual patches
% YrA:     Residual signal at the level of each component
%
% Author: Eftychios A. Pnevmatikakis, Simons Foundation, 2015, 2016

if nargin < 6
    options = CNMFSetParms();
else
    options = CNMFSetParms(options);
end

memmaped = isobject(data);
if ~memmaped
    Y = data;
    clear data;  % TODO check if necessary
    sizY = size(Y);
    Yr = reshape(Y,[],sizY(end));
    F_dark = min(Yr(:));

    % create a memory mapped object named data_file.mat and open it read-only
    if options.create_memmap
        save('data_file.mat','Yr','Y','F_dark','sizY','-v7.3');
        data = matfile('data_file.mat','Writable',false);
        memmaped = true;
    elseif isa(Yr, 'single') || isa(Yr, 'double')
        data = Yr;
    else
        data = single(Yr);
    end
else
    sizY = data.sizY;
    if ismember('F_dark',who(data))
        F_dark = data.F_dark;
    elseif ismember('nY',who(data))
        F_dark = data.nY;
    else
        F_dark = 0;
    end
end
F_dark = double(F_dark);

% ensure correct dimensions options
options.d1 = sizY(1);
options.d2 = sizY(2);
if length(sizY) == 4
    options.d3 = sizY(3);
end

if nargin < 5 || isempty(p)
    p = 0;
end

if nargin < 4 || isempty(tau)
    tau = 5;
end

if nargin < 3 || isempty(patches)
    patches = construct_patches(sizY(1:end-1),[50,50]);  % TODO fix default for 3d case
end
n_patches = length(patches);

if nargin < 2 || isempty(K)
    K = 10;
end

%% running CNMF on each patch, in parallel
%RESULTS(n_patches) = struct('A', [], 'b', [], 'C', [], 'f', [], 'S', [], 'P', []);
RESULTS = CNMF();
RESULTS(n_patches) = CNMF();

if memmaped    
    parfor i = 1:n_patches
        patch_idx = patch_to_indices(patches{i});
        Yp = data.Y(patch_idx{:},:);
        RESULTS(i) = process_patch_object(Yp,F_dark, K, p, tau, options);
        fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
        RESULTS(i).Y = [];
        RESULTS(i).Yr = [];
    end

else  % avoid copying the entire dataset to each worker, for in-memory data
    for i = n_patches:-1:1
        patch_idx = patch_to_indices(patches{i});
        Yp = Y(patch_idx{:},:);
        future_results(i) = parfeval(@process_patch_object, 1, Yp, F_dark, K, p, tau, options);
    end
    for i = 1:n_patches
        [idx, value] = fetchNext(future_results);
        RESULTS(idx) = value;
        RESULTS(idx).Y = [];
        RESULTS(idx).Yr = [];
        fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
    end
end

%% combine results into one structure
fprintf('Combining results from different patches... \n');
P.sn = zeros(sizY(1:end-1));
P.b = {};
P.c1 = {};
P.gn = {};
P.neuron_sn = {};
if isfield(RESULTS(1).P,'sn_ds')
    P.sn_ds = zeros(sizY(1:end-1));
end

cnt = 0;
d = prod(sizY(1:end-1));
A = sparse(d,n_patches*K);
B = sparse(d,n_patches);
MASK = zeros(sizY(1:end-1));
F = zeros(n_patches,sizY(end));
for i = 1:n_patches
    patch_lin_idx = patch_to_linear(patches{i}, sizY);
    patch_size = patches{i}(2:2:end) - patches{i}(1:2:end) + 1;
    for k = 1:K
        if k > size(RESULTS(i).A,2)
            break;
        end
        cnt = cnt + 1;
        A(patch_lin_idx,cnt) = RESULTS(i).A(:,k);
    end
    B(patch_lin_idx,i) = RESULTS(i).b;
    MASK(patch_lin_idx) = MASK(patch_lin_idx) + 1;
    P.sn(patch_lin_idx) = reshape(RESULTS(i).P.sn,patch_size);
    if isfield(RESULTS(i).P,'sn_ds')
        P.sn_ds(patch_lin_idx) = reshape(RESULTS(i).P.sn_ds,patch_size);
    end
    P.b = [P.b;RESULTS(i).P.b];
    P.c1 = [P.c1;RESULTS(i).P.c1];
    P.gn = [P.gn;RESULTS(i).P.gn];
    P.neuron_sn = [P.neuron_sn;RESULTS(i).P.neuron_sn];
    F(i,:) = RESULTS(i).f;
end
A(:,cnt+1:end) = [];
A = spdiags(1./MASK(:),0,d,d)*A;
B = spdiags(1./MASK(:),0,d,d)*B;
C = cell2mat({RESULTS(:).C}');
S = cell2mat({RESULTS(:).S}');
ff = find(sum(A,1)==0);
A(:,ff) = [];
C(ff,:) = [];
S(ff,:) = [];
fprintf(' done. \n');

% estimate active pixels
if options.cluster_pixels
    fprintf('Classifying pixels...')
    P.active_pixels = zeros(sizY(1:end-1));
    psdx_size = [patches{end}(2:2:end), size(RESULTS(1).P.psdx,2)];
    P.psdx = zeros(psdx_size);

    for i = 1:n_patches
        patch_idx = patch_to_indices(patches{i});
        patch_size = patches{i}(2:2:end) - patches{i}(1:2:end) + 1;
        P.active_pixels(patch_idx{:}) = P.active_pixels(patch_idx{:}) + reshape(RESULTS(i).P.active_pixels,patch_size);
        P.psdx(patch_idx{:},:) = reshape(RESULTS(i).P.psdx,[patch_size, psdx_size(end)]);
    end

    if length(sizY) == 3
        X = P.psdx(:,:,1:min(size(P.psdx,3),500));
    else
        X = P.psdx(:,:,:,1:min(size(P.psdx,4),500));
    end
    X = reshape(X,[],size(X,ndims(X)));
    X = bsxfun(@minus,X,mean(X,2));  % center
    X = spdiags(std(X,[],2)+1e-5,0,size(X,1),size(X,1))\X;
    [L,Cx] = kmeans_pp(X',2);
    [~,ind] = min(sum(Cx(max(1,end-49):end,:),1));
    P.active_pixels = (L==ind);
    P.centroids = Cx;
    fprintf(' done. \n');
end

%% merge results
fprintf('Merging overlaping components...')
Am = A;
Cm = C;
Pm = P;
Sm = S;
Km = 0;
Kn = size(A,2);

while Km < Kn
    Kn = size(Am,2);
    [Am,Cm,~,~,Pm,Sm] = merge_components([],Am,[],Cm,[],Pm,Sm,options);
    Km = size(Am,2);
end
fprintf(' done. \n');

%% compute spatial and temporal background using a rank-1 fit
fprintf('Computing background components...')
% fin = [mean(F,1);rand(options.gnb-1,size(F,2))];
% for iter = 1:150
%     fin = diag(sqrt(sum(fin.^2,2)))\fin;
%     bin = full(max(B*(F*fin')/(fin*fin'),0));
%     fin = max((bin'*bin)\(bin'*B)*F,0);
% end
[bin,fin] = fast_nmf(B,F,options.gnb,100);

fprintf(' done. \n');

%% classify components
options.classify_comp = false; % components are now classified within each patch
if options.classify_comp
    fprintf('Classifying components...')
    options.space_thresh = options.patch_space_thresh;
    options.time_thresh = options.patch_time_thresh;
    if ~memmaped
        [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(Y,Am,Cm,bin,fin,options);
    else
        [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(data,Am,Cm,bin,fin,options);
    end
    ind = ind_space & ind_time;
    fprintf(' done. \n');
else
    ind = true(size(Am,2),1);
    rval_space = NaN(size(Am,2),1);
    rval_time = NaN(size(Am,2),1);
end

A = Am(:,ind);
C = Cm(ind,:);

Pm.rval_space = rval_space;
Pm.rval_time = rval_time;
Pm.A_throw = Am(:,~ind);
Pm.C_throw = Cm(~ind,:);

%% update again before screening
if options.refine_flag
    % update spatial components
    fprintf('Updating spatial components...');
    options.nb = options.gnb;
    if ~isfield(Pm,'mis_values'); Pm.mis_values = []; end
    if ~isfield(Pm,'mis_entries'); Pm.mis_entries = []; end
    [A,b,C,Pm] = update_spatial_components(data,C,fin,[A,bin],Pm,options);
    fprintf(' done. \n');

    % update temporal components
    fprintf('Updating temporal components... ')
    Pm.p = 0;
    [C,f,P,S,YrA] = update_temporal_components_fast(data,A,b,C,fin,Pm,options);
    fprintf(' done. \n');
else
    b = bin;
    f = fin;
    AY = mm_fun([A,double(b)],data);
    AA = [A,double(b)]'*[A,double(b)];
    YrA = bsxfun(@times, 1./sum([A,double(b)].^2)',AY - AA*[C;f]);
    YrA = YrA(1:size(C,1),:);
end
P.p = p;
end

function idx = patch_to_indices(patch)
    % helper function to build indices vector from patch start/stop indices
    idx = arrayfun(@(x,y) x:y, patch(1:2:end), patch(2:2:end), 'un', false);
end

function idx = patch_to_linear(patch, sizY)
    % helper function to build linear indices from patch start/stop indices
    slice_idx = patch_to_indices(patch);
    subs_idx = cell(1, numel(slice_idx));
    [subs_idx{:}] = ndgrid(slice_idx{:});
    subs_idx = cellfun(@(x) x(:), subs_idx, 'un', false);
    idx = sub2ind(sizY(1:end-1), subs_idx{:});
end

function result = process_patch(Y, F_dark, K, p, tau, options)
    % helper function to apply CNMF to a small patch

    sizY = size(Y);
    options.d1 = sizY(1);
    options.d2 = sizY(2);
    if ndims(Y) == 3
        options.d3 = 1;
    else
        options.d3 = sizY(3);
    end
    options.nb = 1;
    options.temporal_parallel = 0;  % turn off parallel updating for temporal components
    options.spatial_parallel = 0;   % turn off parallel updating for spatial components
    options.space_thresh = options.patch_space_thresh;    % put a low acceptance threshold initially
    options.time_thresh = options.patch_time_thresh;

    Y = double(Y - F_dark);
    Y(isnan(Y)) = F_dark;

    [P,Y] = preprocess_data(Y,p,options);
    Yr = reshape(Y,[],sizY(end));

    [Ain,Cin,bin,fin] = initialize_components(Y,K,tau,options,P);
    [A,b,Cin,P] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
    P.p = 0;
    options.p = 0;
    [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

    if ~isempty(A) && ~isempty(C)
        [Am,Cm,~,~,P] = merge_components(Yr,A,b,C,f,P,S,options);
        [A,b,Cm,P] = update_spatial_components(Yr,Cm,f,[Am,b],P,options);
        [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cm,f,P,options);
        [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(Y,A,C,b,f,options); 
        %ind = ind_space & ind_time;        
        ind_corr = ind_space;
        
        try  % matlab 2017b or later is needed
            [ind_cnn,value] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
        catch
            ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
        end     

        fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
        ind_exc = (fitness < options.min_fitness);
        ind = (ind_corr | ind_cnn) & ind_exc;
        %fitness_delta = compute_event_exceptionality(diff(C+YrA,[],2),0);
        %ind = (ind_space & ind_time) | (fitness < options.patch_max_fit) | (fitness_delta < options.patch_max_fit_delta);
        
        P.rval_space = rval_space;
        P.rval_time = rval_time;
        P.ind_space = ind_space;
        P.ind_time = ind_time;
        P.fitness = fitness;
        P.fitness_delta = fitness_delta;
        P.A_throw = A(:,~ind);
        P.C_throw = C(~ind,:);
    end
    result.A = A(:,ind);
    result.b = b;
    result.C = C(ind,:);
    result.f = f;
    result.S = S;
    result.P = P;
end

function CNM = process_patch_object(Y,F_dark,K,p,tau,options)
    CNM = CNMF();
    if ndims(Y) > 3; d3 = size(Y,3); else; d3 = 1; end
    options = CNMFSetParms(options,...
                'd1',size(Y,1),...
                'd2',size(Y,2),...
                'd3',d3,...
                'p',p,...
                'gSig',tau,...
                'temporal_parallel',false,...
                'spatial_parallel',false,...
                'space_thresh',options.patch_space_thresh,...
                'time_thresh',options.patch_time_thresh,...
                'cnn_thr',options.patch_cnn_thr,...
                'min_fitness',options.patch_min_fitness);
    Y = single(Y) - single(F_dark);
    Y(isnan(Y)) = single(F_dark);
    CNM.fit(Y,options,K);                
end