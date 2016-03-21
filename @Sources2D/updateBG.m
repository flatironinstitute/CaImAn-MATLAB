function [B, F] = updateBG(obj, Y, nb, model)
%% initialize background given the rank of background
% input:
%   Y: d*T or d1*d2*T, imaging data
%   nb: scalar, number of the background
%   mdoel: string, {'nmf', 'svd'}
% output:
%   B:  d*nb matrix, spatial components of the background
%   F:  nb*T matrix, temporal components of the background

%% parameters
if ~exist('nb', 'var')
    nb = obj.options.nb;
else
    obj.options.nb = nb;
end

if ~exist('model', 'var'); model = 'nmf'; end

% convert data into vectorized format, each frame (d1*d2) is one column
Y = obj.reshape(Y, 1);
[~, T] = size(Y,2);

%% start factorization
if strcmpi(model, 'svd')
    Y_mean = mean(Y, 1);
    Y = bsxfun(@minus, Y,Y_mean);
    [~, ~, v] = svdsecon(Y(1:10:end, :), nb);
    B = [ones(d1*d2,1), Y*v];
    F = [Y_mean; v'];
else
    tsub = max(1, round(T/500));   % maximumly 1000 frames are used for initialization
    frames = (1:tsub:T);
    [B, f] = nnmf(Y(:, frames), nb); %run nmf on subsampled data for initialization
    F = imresize(f, [nb, T]);
    for miter=1:5
        % update temporally
        U = B'*Y;
        V = B'*B;
        bb = diag(V);
        for m=1:nb
            F(m, :) = max(0, F(m,:)+(U(m, :)-V(m,:)*F)/bb(m));
        end
        
        % update spatially
        U = Y*F';
        V = F*F';
        ff = diag(V);
        for m=1:nb
            B(:, m) = max(0, B(:, m)+(U(:, m)- B*V)/ff(m));
        end
    end
end
obj.b = B;
obj.f = F;

function [U,S,V] = svdsecon(X,k)
% Input:
% X : m x n matrix
% k : extracts the first k singular values
%
% Output:
% X = U*S*V' approximately (up to k)
%
% Description:
% Does equivalent to svds(X,k) but faster
% Requires that k < min(m,n) where [m,n] = size(X)
% This function is useful if k is much smaller than m and n
% or if X is sparse (see doc eigs)
%
% Vipin Vijayan (2014)

%X = bsxfun(@minus,X,mean(X,2));
[m,n] = size(X);
assert(k <= m && k <= n, 'k needs to be smaller than size(X,1) and size(X,2)');

if  m <= n
    C = X*X';
    [U,D] = eigs(C,k);
    clear C;
    if nargout > 2
        V = X'*U;
        s = sqrt(abs(diag(D)));
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
    C = X'*X; 
    [V,D] = eigs(C,k);
    clear C;
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(abs(diag(D)));
    U = bsxfun(@(x,c)x./c, U, s');
    S = diag(s);
end
