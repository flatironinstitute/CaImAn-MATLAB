function C = HALS_temporal(Y, A, C, maxIter,tol)
%% run HALS by fixating all spatial components 
% input: 
%   Y:       d*T, fluorescence data
%   A:       d*K, spatial components 
%   C:       K*T, temporal components 
%   maxIter: maximum number of iterations
%   tol:     convergence threshold 
% output: 
%   C: K*T, updated temporal components 

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes
% Friedrich's NIPS paper "Fast Constrained Non-negative Matrix
% Factorization for Whole-Brain Calcium Imaging Data

%% options 
if nargin < 5; tol = 1e-4; end
if nargin < 4 || isempty(maxIter); maxIter=1; end

%% initialization 
K = size(A, 2);     % number of components 
U = A'*Y; 
V = A'*A; 
aa = diag(V);   % squares of l2 norm all all components 
%% updating

repeat = 1;
miter = 0;
while repeat
    C_ = C;
    for k=1:K
        ck = max(0, C(k, :) + (U(k, :)-V(k, :)*C)/aa(k));     % update ck
        C(k, :) = ck; 
    end
    repeat = (miter < maxIter) && norm(C - C_,'fro')/norm(C_,'fro') > tol;
    miter = miter + 1;
end