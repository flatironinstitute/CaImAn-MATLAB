function ff = classify_components(A,P,options)

% Classify components based on their spatial overlap with the active pixels
% Each component is classified as true if when restricted to the active
% pixels retains at least options.cl_thr of its l_2 norm

% INPUTS
% A         : set of spatial components (d x K)
% P         : parameter structure (must include P.active_pixels) obtained from preprocess_data.m
% options   : options structure (optional) for reading the energy threshold

% OUTPUT
% ff        : K x 1 binary vector classifying the components

% Author: Eftychios A. Pnevmatikakis
%           Simons Foundation, 2016

defoptions = CNMFSetParms;
if nargin < 3 || isempty(options)
    options = defoptions;
end

if ~isfield(options,'cl_thr') || isempty(options.cl_thr)
    cl_thr = defoptions.cl_thr;
else
    cl_thr = options.cl_thr;
end

ff = false(1,size(A,2));
for i = 1:size(A,2)
    a1 = A(:,i);
    a2 = A(:,i).*P.active_pixels(:);
    if sum(a2.^2) >= cl_thr^2*sum(a1.^2)
        ff(i) = true;
    end
end