function [A, C, b, f] = HALS(Y, A, C, b, f, params)
%% Hierarchical alternating least square method for solving NMF problem  
% Y = A*C + b*f 

%input: 
%   Y:      d1 X d2 X T, raw data. It will be reshaped to (d1*d2) X T in this
%   function 
%   A:      (d1*d2) X K, initial value of spatial components 
%   C:      K X T, initial value of temporal components 
%   b:      (d1*d2) X 1, initial value of background spatial component 
%   f:      1 X T, initial value of background temporal component
%   params: parameters used in this function. 
%       bSiz:   blur size. A box kernel (bSiz X bSiz) will be convolved
%       with each neuron's initial spatial component, then all nonzero
%       pixels will be picked as pixels to be updated, and the rest will be
%       forced to be 0. 
%       maxIter: maximum iteration of iterating HALS. 

% Author: Pengcheng Zhou, Columbia University, based on a python
% implementation from Johannes Friedrich, Columbia University, 2015.

%[d1, d2, ~] = size(Y);
dimY = ndims(Y)-1;
sizY = size(Y);
d = sizY(1:dimY);

if nargin < 6; params = []; end
if isfield(params, 'bSiz'), bSiz = params.bSiz; else bSiz=3; end
if isfield(params, 'maxIter'), maxIter = params.maxIter; else maxIter=5; end

blur_kernel = ones(bSiz);
if dimY == 3; blur_kernel(:,:,3:end) = []; end
blur_kernel = blur_kernel/numel(blur_kernel);
if ismatrix(A)
    A = reshape(A,[d,size(A,2)]);
    ind_A = (imfilter(A,blur_kernel)>0);
end
A = sparse(reshape(A,prod(d),[]));
ind_A = sparse(reshape(ind_A,prod(d),[]));
K = size(A,2);
% if ndims(A)==3
%     ind_A = (imfilter(A, blur_kernel)>0); 
%     A = reshape(A, d1*d2, []);
%     ind_A = reshape(ind_A, d1*d2, []); 
% else 
%     ind_A = (imfilter(reshape(A, d1,d2,[]), blur_kernel)>0); 
%     ind_A = reshape(ind_A, d1*d2, []); 
% end
% ind_A = sparse(ind_A);  %indicator of nonnero pixels 
% K = size(A, 2);         %number of neurons 

%% update spatial and temporal components neuron by neurons
%Yres = reshape(Y, d1*d2, []) - A*C - b*f;
Yres = reshape(Y, prod(d), []) - A*C - b*f;

for miter=1:maxIter
    for mcell = 1:K
        ind_pixels = ind_A(:, mcell);
        tmp_Yres = Yres(ind_pixels, :);
        
        % update temporal component of the neuron
        c0 = C(mcell, :);
        a0 = A(ind_pixels, mcell);
        norm_a2 = norm(a0, 2)^2;
        C(mcell, :) = max(0, c0 + a0'*tmp_Yres/norm_a2);
        tmp_Yres = tmp_Yres + a0*(c0-C(mcell,:));
        
        % update spatial component of the neuron
        norm_c2 = norm(C(mcell,:),2)^2;
        tmp_a = max(0, a0 + tmp_Yres*C(mcell, :)'/norm_c2);
        A(ind_pixels, mcell) = tmp_a; 
        Yres(ind_pixels,:) = tmp_Yres + (a0-tmp_a)*C(mcell,:);
    end
    
    % update temporal component of the background
    f0 = f;
    b0 = b;
    norm_b2 = norm(b0,2)^2;
    f = max(0, f0 + b0'*Yres/norm_b2);
    Yres = Yres + b0*(f0-f);
    % update spatial component of the background
    norm_f2 = norm(f, 2)^2;
    b = max(0, b0 + Yres*f'/norm_f2);
    Yres = Yres + (b0-b)*f;
    
    %fprintf('Iteration %d, the norm of  residual is %.4f\n', miter, norm(Yres, 'fro'));
end