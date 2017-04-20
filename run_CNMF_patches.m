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

    % create a read-only memory mapped object named data_file.mat
    if options.create_memmap
        save('data_file.mat','Yr','Y','F_dark','sizY','-v7.3');
        data = matfile('data_file.mat','Writable',false);
        memmaped = true;
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

Yc = cell(n_patches,1);
if ~memmaped
    for i = 1:n_patches
        patch_idx = patch_to_indices(patches{i});
        Yc{i} = Y(patch_idx{:}, :);
    end
end

if nargin < 2 || isempty(K)
    K = 10;
end

%% running CNMF on each patch, in parallel
RESULTS(n_patches) = struct('A', [], 'b', [], 'C', [], 'f', [], 'S', [], 'P', []);
parfor i = 1:n_patches
    if memmaped
        patch_idx = patch_to_indices(patches{i});
        Yp = data.Y(patch_idx{:},:);
    else
        Yp = Yc{i};
    end
    if length(sizY) == 3
        [d1,d2,~] = size(Yp);
        d3 = 1;
    else
        [d1,d2,d3,~] = size(Yp);
    end

    RESULTS(i) = process_patch(Yp, [d1, d2, d3], F_dark, K, p, tau, options);
    fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
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

if options.cluster_pixels
    P.active_pixels = zeros(sizY(1:end-1));
    if length(sizY) == 3
        P.psdx = zeros(patches{end}(2),patches{end}(4),size(RESULTS(1).P.psdx,2));
    else
        P.psdx = zeros(patches{end}(2),patches{end}(4),patches{end}(6),size(RESULTS(1).P.psdx,2));
    end
end

cnt = 0;
d = prod(sizY(1:end-1));
A = sparse(d,n_patches*K);
B = sparse(d,n_patches);
MASK = zeros(sizY(1:end-1));
F = zeros(n_patches,sizY(end));
for i = 1:n_patches
    patch_idx = patch_to_indices(patches{i});
    patch_size = patches{i}(2:2:end) - patches{i}(1:2:end) + 1;
    for k = 1:K
        if k > size(RESULTS(i).A,2)
            continue
        end
        cnt = cnt + 1;
        Atemp = zeros(sizY(1:end-1));
        Atemp(patch_idx{:}) = reshape(full(RESULTS(i).A(:,k)),patch_size);
        A(:,cnt) = sparse(Atemp(:));
    end
    b_temp = zeros(sizY(1:end-1));
    b_temp(patch_idx{:}) = reshape(full(RESULTS(i).b),patch_size);
    MASK(patch_idx{:}) = MASK(patch_idx{:}) + 1;
    P.sn(patch_idx{:}) = reshape(RESULTS(i).P.sn,patch_size);
    if isfield(RESULTS(i).P,'sn_ds')
        P.sn_ds(patch_idx{:}) = reshape(RESULTS(i).P.sn_ds,patch_size);
    end
    if options.cluster_pixels
        P.active_pixels(patch_idx{:}) = P.active_pixels(patch_idx{:}) + reshape(RESULTS(i).P.active_pixels,patch_size);
        patch_size_cell = num2cell(patch_size);
        P.psdx(patch_idx{:},:) = reshape(RESULTS(i).P.psdx,patch_size_cell{:},[]);
    end
    P.b = [P.b;RESULTS(i).P.b];
    P.c1 = [P.c1;RESULTS(i).P.c1];
    P.gn = [P.gn;RESULTS(i).P.gn];
    P.neuron_sn = [P.neuron_sn;RESULTS(i).P.neuron_sn];
    B(:,i) = sparse(b_temp(:));
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

if options.cluster_pixels
    % estimate active pixels
    fprintf('Classifying pixels...')
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
fin = [mean(F);rand(options.gnb-1,length(F))];
for iter = 1:150
    fin = diag(sqrt(sum(fin.^2,2)))\fin;
    bin = max(B*(F*fin')/(fin*fin'),0);
    fin = max((bin'*bin)\(bin'*B)*F,0);
end
fprintf(' done. \n');

%% classify components
if options.classify_comp
    fprintf('Classifying components...')
    options.space_thresh = 0.3;
    options.time_thresh = 0.3;
    options.max_pr_thr = 0.75;
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

%% update spatial components
fprintf('Updating spatial components...');
options.nb = options.gnb;
if ~isfield(Pm,'mis_values'); Pm.mis_values = []; end
if ~isfield(Pm,'mis_entries'); Pm.mis_entries = []; end
[A,b,C,Pm] = update_spatial_components(data,C,fin,[A,bin],Pm,options);
fprintf(' done. \n');

%% update temporal components
fprintf('Updating temporal components... ')
Pm.p = 0;
[C,f,P,S,YrA] = update_temporal_components_fast(data,A,b,C,fin,Pm,options);
fprintf(' done. \n');

end

function idx = patch_to_indices(patch)
    % helper function to build indices vector from patch start/stop indices
    idx = arrayfun(@(x,y) x:y, patch(1:2:end), patch(2:2:end), 'un', false);
end

function result = process_patch(Y, dims, F_dark, K, p, tau, options)
    % helper function to apply CNMF to a small patch

    options.d1 = dims(1);
    options.d2 = dims(2);
    options.d3 = dims(3);
    options.nb = 1;
    options.temporal_parallel = 0; % turn off parallel updating for temporal components
    options.spatial_parallel = 0;  % turn off parallel updating for spatial components

    Y = double(Y - F_dark);
    Y(isnan(Y)) = F_dark;

    [P,Y] = preprocess_data(Y,p);
    Yr = reshape(Y,prod(dims),[]);

    [Ain,Cin,bin,fin] = initialize_components(Y,K,tau,options,P);
    [A,b,Cin,P] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
    P.p = 0;
    [C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

    if ~isempty(A) && ~isempty(C)
        [Am,Cm,~,~,P] = merge_components(Yr,A,b,C,f,P,S,options);
        [A,b,Cm,P] = update_spatial_components(Yr,Cm,f,[Am,b],P,options);
        P.p = p;
        [C,f,P,S] = update_temporal_components(Yr,A,b,Cm,f,P,options);
    end

    result.A = A;
    result.b = b;
    result.C = C;
    result.f = f;
    result.S = S;
    result.P = P;
end