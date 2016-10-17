function A = nnls_spatial_thresh(Y, A, C, active_pixel, maxN, smin, sn)
%% run HALS by fixating all spatial components
% input:
%   Y:  d*T, fluorescence data
%   A:  d*K, spatial components
%   C:  K*T, temporal components
%   active_pixel: d*T binary matrix, indicating the nonzero elements of A
%   maxN: scalar, maximum number of neurons overlapping at one pixel
% output:
%   A: d*K, updated spatial components

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes

%% options for HALS
if nargin<5;    maxN = 5;    end;    %maximum iteration number
if nargin<4;    active_pixel=true(size(A));
elseif isempty(active_pixel)
    active_pixel = true(size(A));
else
    active_pixel = logical(active_pixel);
end;     %determine nonzero pixels
ind_fit = find(sum(active_pixel,2)>1e-9);


%% initialization
Ymean = mean(Y,2); 
Y(bsxfun(@lt, Y, Ymean+smin*sn)) = 0; 
C(C<smin) = 0; 
CC = C*C';
YC = C*Y';
A = zeros(size(A)); 

%% updating
for m=1:length(ind_fit)
    ind = active_pixel(ind_fit(m), :);
    A(ind_fit(m), ind) = nnls(CC(ind, ind), YC(ind, ind_fit(m)), [], 1e-4, maxN);
end
function s = nnls(A, b, s, tol, maxIter)
%% fast algorithm for solving nonnegativity constrained least squared
% problem minize norm(y-K*s, 2), s.t. s>=0.

%% inputs:
%   A: n x p matrix, K'*K
%   b: n x 1 vector, K'*y
%   s: p x 1 vector, warm started s
%   tol: scalar, smallest nonzero values
%   maxIter: scalar, maximum nonzero values

%% outputs:
%   s: p x 1 vector, solution

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging
% Bro R & Jong S, Journal of Chemometrics 1997, A FAST NON-NEGATIVITY-CONSTRAINED LEAST SQUARES ALGORITHM

%% input arguments
p = size(A,2);      % dimension of s
if ~exist('s', 'var') || isempty(s)
    s = zeros(p, 1);
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-9;
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = p;
end
if sum(s>0)>maxIter
    s = zeros(p,1);
end
for miter=1:maxIter
    l = b - A*s;            % negative gradient
    Pset = (s>0);       % passive set
    
    if max(l) < tol         % no more passive set
        break;
    end
    
    [~, temp] = max(l);         % choose the one with the largest gradient
    Pset(temp) = true;         % add it to the passive set
    if sum(Pset)>maxIter
        break;
    end
    % correct nonnegativity violations
    while any(Pset)
        % run unconstrained least squares for variables in passive sets
        try
            mu = A(Pset, Pset) \ b(Pset);
        catch
            mu = (A(Pset, Pset) + tol*eye(sum(Pset))) \ b(Pset);
        end
        
        if all(mu>tol)
            break;
        end
        
        temp = s(Pset) ./ (s(Pset)-mu);
        temp(mu>tol) = [];
        a = min(temp);
        s(Pset) = s(Pset) + a*(mu-s(Pset));
        Pset(s<tol) = false;
    end
    
    s(Pset) = mu;
end