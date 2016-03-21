function [C, S] = HALS_temporal_constrain(Y, A, C, IND, maxIter, P, options)
%% run HALS by fixating all spatial components
% input:
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components
%   C:  K*T, temporal components
%   maxIter:    maximum iteration number
% output:
%   C: K*T, updated temporal components

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes
% Friedrich's NIPS paper "Fast Constrained Non-negative Matrix
% Factorization for Whole-Brain Calcium Imaging Data

%% options
T = size(Y, 2);
K = size(A, 2);     % number of components
if ~exist('maxIter', 'var');    maxIter=1; end
if ~exist('IND', 'var') || isempty(IND)
    IND = true(size(A));
elseif isscalar(IND)
    IND = bsxfun(@gt, A, max(A, [], 1)*IND);
end

if exist('options', 'var') && (~isempty(options))
    method = options.deconv_method;
    restimate_g = options.restimate_g;
    
    if ~exist('P', 'var') || isempty(P)
        active_pixels = find(sum(A,2));                                 % pixels where the greedy method found activity
        unsaturated_pixels = find_unsaturatedPixels(Y);                 % pixels that do not exhibit saturation
        options.pixels = intersect(active_pixels,unsaturated_pixels);   % base estimates only on unsaturated, active pixels
        P.b = cell(K,1);
        P.c1 = cell(K,1);
        P.neuron_sn = cell(K,1);
    end
else
    method = [];
end
%% initialization
U = (A.*IND)'*Y;
V = (A.*IND)'*A;
aa = sum((A.*IND).^2, 1);   % squares of l2 norm all all components
S = zeros(size(C));

%% updating
for miter=1:maxIter
    for k=1:K
        %% deconvolution
        ck = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);
        
        if strcmpi(method,  'constrained_foopsi')
            if restimate_g
                [cc,cb,c1,gn,~,spk] = constrained_foopsi(ck,[],[],[],[],options);
                P.gn{k} = gn;
            else
                [cc,cb,c1,gn,~,spk] = constrained_foopsi(ck,[],[],P.g,[],options);
            end
            gd = max(roots([1,-gn']));  % decay time constant for initial concentration
            gd_vec = gd.^((0:T-1));
            C(k,:) = full(cc(:)' + cb + c1*gd_vec);
            S(k,:) = spk(:)';
        else
            C(k, :) = max(0, ck);
        end
    end
end

