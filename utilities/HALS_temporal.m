function C = HALS_temporal(Y, A, C, maxIter)
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
if nargin<4;    maxIter=1; end

%% initialization 
K = size(A, 2);     % number of components 
U = A'*Y; 
V = A'*A; 
aa = diag(V);   % squares of l2 norm all all components 
%% updating 
for miter=1:maxIter
    for k=1:K
        ck = max(0, C(k, :) + (U(k, :)-V(k, :)*C)/aa(k));     % update ck
        C(k, :) = ck; 
    end
end

