function [A,C,b,f] = HALS_initialization(Y,K,options)

defoptions = CNMFSetParms;
if ~isfield(options,'max_iter_hals_in') || isempty(options.max_iter_hals_in); options.max_iter_hals_in = defoptions.max_iter_hals_in; end
if ~isfield(options,'rem_prct') || isempty(options.rem_prct); options.rem_prct = defoptions.rem_prct; end
if ~isfield(options,'nb') || isempty(options.nb); nb = defoptions.nb; else nb = options.nb; end

init_hals_method = 'random';

if ~exist('K', 'var'),  %find K neurons
    K = 30;  
    warning(['number of neurons are not specified, set to be the default value', num2str(K)]);
end

dimY = ndims(Y) - 1;  % dimensionality of imaged data (2d or 3d)
sizY = size(Y);
T = sizY(end);        % # of timesteps
dx = sizY(1:dimY);    % # of voxels in each axis
d = prod(dx);         % total # of voxels  

med = prctile(Y,options.rem_prct,ndims(Y));
Y = bsxfun(@minus, Y, med);
if strcmpi(init_hals_method,'random');
    A = rand(d,K);
    Y = reshape(Y,d,T);
    C = HALS_temporal(Y,A,rand(K,T),10);
elseif strcmpi(init_hals_method,'cor_im');
    sk = max(round(T/1000),1);
    Cn = correlation_image(Y(:,:,1:sk:T));
    Y = reshape(Y,d,T);
    Cnf = imgaussfilt(Cn,2.25);
    BW = imregionalmax(Cnf);
    C = Y(BW(:),:);
    A = max(Y*pinv(C),0);
end

for iter = 1:options.max_iter_hals_in
    A = HALS_spatial(Y, A, C);
    C = HALS_temporal(Y, A, C);
end

ind_del = find(std(C,0,2)==0); 
A(:, ind_del) = []; 
C(ind_del, :) = []; 

Y = bsxfun(@plus,Y-A*C,med(:));

b = [med(:),rand(d,nb-1)];

for nmfiter = 1:100
    f = max((b'*b)\(b'*Y),0);
    b = max((Y*f')/(f*f'),0);    
end