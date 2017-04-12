function A = HALS_spatial(Y, A, C, active_pixel, maxIter)
%% run HALS by fixating all spatial components 
% input: 
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components 
%   C:  K*T, temporal components 
% output: 
%   A: d*K, updated spatial components 

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes
% Friedrich's NIPS paper "Fast Constrained Non-negative Matrix
% Factorization for Whole-Brain Calcium Imaging Data

%% options for HALS
if nargin<5;    maxIter = 1;    end;    %maximum iteration number 
if nargin<4;    active_pixel=true(size(A));
elseif isempty(active_pixel)
    active_pixel = true(size(A)); 
end;     %determine nonzero pixels 

%% initialization 
A(~active_pixel) = 0; 
K = size(A, 2);     % number of components 
U = mm_fun(C,Y); 
V = C*C'; 
cc = diag(V);   % squares of l2 norm all all components 

%% updating 
for miter=1:maxIter
    for k=1:K
        tmp_ind = active_pixel(:, k); 
        if any(tmp_ind)
            ak = max(0, A(tmp_ind, k)+(U(tmp_ind, k)-A(tmp_ind,:)*V(:, k))/cc(k)); 
            A(tmp_ind, k) = ak; 
        end
    end
end