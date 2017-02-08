function [C,U] = HALS_temporal(Y, A, C, maxIter,tol,nonneg,flag_AY)
%% run HALS by fixating all spatial components 
% input: 
%   Y:       d*T, fluorescence data or K*T, product A'*Y if flag_AY = true 
%   A:       d*K, spatial components 
%   C:       K*T, temporal components 
%   maxIter: maximum number of iterations
%   tol:     convergence threshold 
%   noneg:   restrict output to be non-negative
% output: 
%   C: K*T, updated temporal components 
%   U:      product A'*Y (for downstream use)

% Author: Pengcheng Zhou, Carnegie Mellon University and Eftychios A. Pnevmatikakis, 
% adapted from Johannes Friedrich's NIPS paper "Fast Constrained Non-negative Matrix
% Factorization for Whole-Brain Calcium Imaging Data

%% options 
if nargin < 7 || isempty(flag_AY); flag_AY = false; end
if nargin < 6 || isempty(nonneg); nonneg = true; end
if nargin < 5 || isempty(tol); tol = 1e-4; end
if nargin < 4 || isempty(maxIter); maxIter=1; end

%% initialization 
K = size(A, 2);     % number of components 
if flag_AY
    U = Y;
else
    U = mm_fun(A,Y);
end
 
V = A'*A; 
aa = diag(V);   % squares of l2 norm all all components 
%% updating

repeat = 1;
miter = 0;
while repeat
    C_ = C;
    for k=1:K
        ck = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);
        if nonneg
            ck = max(0, ck);     % update ck
        end
        C(k, :) = ck; 
    end
    repeat = (miter < maxIter) && norm(C - C_,'fro')/norm(C_,'fro') > tol;
    miter = miter + 1;
end