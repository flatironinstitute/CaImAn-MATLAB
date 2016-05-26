function  merged_ROIs = quickMerge(obj, temporal_component)
%% merge neurons based on simple spatial and temporal correlation 
% input: 
%   temporal_component:     character, {'C', 'S'}, it compupte temporal
%   correlation based on either calcium traces ('C') or spike counts ('S').
% output: 
%   merged_ROIs: cell arrarys, each element contains indices of merged
%   components 

%% Author: Pengcheng Zhou, Carnegie Mellon University. 
%  The basic idea is proposed by Eftychios A. Pnevmatikakis: high temporal
%  correlation + spatial overlap
%  reference: Pnevmatikakis et.al.(2016). Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data. Neuron

%% variables & parameters
A = obj.A;          % spatial components 
C = obj.C;          % temporal components 
options = obj.options;      % options for merging 
merge_thr = options.merge_thr;      % merging threshold 
[K, ~] = size(C);   % number of neurons 

%% find neuron pairs to merge
% compute spatial correlation
A_corr = corr(A);

% compute temporal correlation
if ~exist('temporal_component', 'var') && strcmpi(temporal_component, 'S')
    C_corr = corr(obj.S) - eye(K); 
else
    C_corr = corr(C')-eye(K);
end

%% using merging criterion to detect paired neurons
flag_merge = and((C_corr>=merge_thr), A_corr>0.1);

[l,c] = graph_connected_comp(sparse(flag_merge));     % extract connected components

MC = bsxfun(@eq, reshape(l, [],1), 1:c);
MC(:, sum(MC,1)==1) = [];
if isempty(MC)
    fprintf('All pairs of neurons are below the merging criterion!\n\n');
    merged_ROIs = []; 
    return;
else
    fprintf('%d neurons will be merged into %d new neurons\n\n', sum(MC(:)), size(MC,2));
end

%% start merging
[nr, n2merge] = size(MC);
ind_del = false(nr, 1 );    % indicator of deleting corresponding neurons
merged_ROIs = cell(n2merge,1); 

for m=1:n2merge
    IDs = find(MC(:, m));   % IDs of neurons within this cluster
    merged_ROIs{m} = IDs; 
    
    % determine searching area
    active_pixel = (sum(A(:,IDs), 2)>0);  
    
    % update spatial/temporal components of the merged neuron
    data = A(active_pixel, IDs)*C(IDs, :); 
    ci = C(IDs(1), :); 
    for miter=1:10
        ai = data*ci'/(ci*ci'); 
        ci = ai'*data/(ai'*ai); 
    end 
    % normalize ai to make its maximum to be 1
    max_ai = max(ai, [], 1);           
    A(active_pixel, IDs(1)) = ai/max_ai;
    C(IDs(1), :) = ci*max_ai; 
    
    % remove merged elements
    ind_del(IDs(2:end)) = true;
end

% remove merged neurons and update obj
A(:, ind_del) = [];
C(ind_del, :) = [];
obj.A = A; 
obj.C = C; 
