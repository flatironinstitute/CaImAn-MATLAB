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
    Yr = reshape(Y,prod(sizY(1:end-1)),[]);
    F_dark = min(Yr(:));

    % create a read-only memory mapped object named data_file.mat
    if options.create_memmap
        save('data_file.mat','Yr','Y','F_dark','sizY','-v7.3');
        data = matfile('data_file.mat','Writable',false);
        memmaped = true;
    else
        data = Yr;
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

if isempty(options.d1); options.d1 = sizY(1); end
if isempty(options.d2); options.d2 = sizY(2); end
if isempty(options.d3)
    if length(sizY) == 3
        options.d3 = 1;
    else
        options.d3 = sizY(3);
    end
end

if nargin < 5 || isempty(p)
    p = 0;
end

if nargin < 4 || isempty(tau)
    tau = 5;
end

if nargin < 3 || isempty(patches)
    patches = construct_patches(sizY(1:end-1),[50,50]);
end

Yc = cell(length(patches),1);
if ~memmaped
    for i = 1:length(patches)
        if length(sizY) == 3
            Yc{i} = Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:);
        else
            Yc{i} = Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6),:);
        end
    end
end        

if nargin < 2 || isempty(K)
    K = 10;
end

RESULTS(length(patches)) = ...
    struct('A', [], 'b', [], 'C', [], 'f', [], 'S', [], 'P', []);

%% running CNMF on each patch, in parallel
n_patches = length(patches);
parfor i = 1:n_patches
    if length(sizY) == 3
        if memmaped
            Yp = data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:);
        else
            Yp = Yc{i};
        end
        [d1,d2,~] = size(Yp);
        d3 = 1;
    else
        if memmaped
            Yp = data.Y(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6),:);
        else
            Yp = Yc{i};
        end
        [d1,d2,d3,~] = size(Yp);
    end

    RESULTS(i) = process_patch(Yp, [d1, d2, d3], F_dark, K, p, tau, options);
    fprintf(['Finished processing patch # ',num2str(i),' out of ',num2str(n_patches), '.\n']);
end

%% combine results into one structure
fprintf('Combining results from different patches... \n');
d = prod(sizY(1:end-1));
A = sparse(d,length(patches)*K);
P.sn = zeros(sizY(1:end-1));
if isfield(RESULTS(1).P,'sn_ds'); P.sn_ds = zeros(sizY(1:end-1)); end
IND = zeros(sizY(1:end-1));
P.b = {};
P.c1 = {};
P.gn = {};
P.neuron_sn = {};

if options.cluster_pixels
    P.active_pixels = zeros(sizY(1:end-1));
    if length(sizY) == 3
        P.psdx = zeros(patches{end}(2),patches{end}(4),size(RESULTS(1).P.psdx,2));
    else
        P.psdx = zeros(patches{end}(2),patches{end}(4),patches{end}(6),size(RESULTS(1).P.psdx,2));
    end
end
    
cnt = 0;
B = sparse(prod(sizY(1:end-1)),length(patches));
MASK = zeros(sizY(1:end-1));
F = zeros(length(patches),sizY(end));
for i = 1:length(patches)
    for k = 1:K
        if k <= size(RESULTS(i).A,2)
            cnt = cnt + 1;
            Atemp = zeros(sizY(1:end-1));
            if length(sizY) == 3
                Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).A(:,k),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);            
            else
                Atemp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = reshape(full(RESULTS(i).A(:,k)),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
            end
            A(:,cnt) = sparse(Atemp(:));
        end
    end
    if length(sizY) == 3
        b_temp = sparse(sizY(1),sizY(2));
        b_temp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).b,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
        P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
        if isfield(RESULTS(i).P,'sn_ds'); P.sn_ds(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = reshape(RESULTS(i).P.sn_ds,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1); end        
        IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + 1;
        if options.cluster_pixels
            P.active_pixels(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) = P.active_pixels(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4)) + ...
                reshape(RESULTS(i).P.active_pixels,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1);
            P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),:) = reshape(RESULTS(i).P.psdx,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,[]);
        end
    else
        b_temp = zeros(sizY(1),sizY(2),sizY(3));
        b_temp(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = reshape(full(RESULTS(i).b),patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
        MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = MASK(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) + 1;
        P.sn(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = reshape(RESULTS(i).P.sn,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
        if isfield(RESULTS(i).P,'sn_ds'); P.sn_ds(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = reshape(RESULTS(i).P.sn_ds,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);  end             
        IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = IND(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) + 1;
        if options.cluster_pixels
            P.active_pixels(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) = P.active_pixels(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6)) + ...
                reshape(RESULTS(i).P.active_pixels,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1,patches{i}(6)-patches{i}(5)+1);
            P.psdx(patches{i}(1):patches{i}(2),patches{i}(3):patches{i}(4),patches{i}(5):patches{i}(6), :) = reshape(RESULTS(i).P.psdx,patches{i}(2)-patches{i}(1)+1,patches{i}(4)-patches{i}(3)+1, patches{i}(6)-patches{i}(5)+1,[]);
        end
    end
    P.b = [P.b;RESULTS(i).P.b];
    P.c1 = [P.c1;RESULTS(i).P.c1];
    P.gn = [P.gn;RESULTS(i).P.gn];
    P.neuron_sn = [P.neuron_sn;RESULTS(i).P.neuron_sn];
    B(:,i) = sparse(b_temp(:));
    F(i,:) = RESULTS(i).f;
end
A(:,cnt+1:end) = [];
A = spdiags(1./MASK(:),0,prod(sizY(1:end-1)),prod(sizY(1:end-1)))*A;
B = spdiags(1./MASK(:),0,prod(sizY(1:end-1)),prod(sizY(1:end-1)))*B;
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
    [rval_space,rval_time,ind_space,ind_time] = classify_comp_corr(data,Am,Cm,bin,fin,options);
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
options.d1 = sizY(1);
options.d2 = sizY(2);
if length(sizY) == 4; options.d3 = sizY(3); end
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