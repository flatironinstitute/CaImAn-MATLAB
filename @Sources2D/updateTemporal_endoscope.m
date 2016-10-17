function [C_offset] = updateTemporal_endoscope(obj, Y, smin)
%% run HALS by fixating all spatial components
% input:
%   Y:  d*T, fluorescence data
%   smin: scalar, threshold for detecting one spikes (>smin*sigma)
% output:
%   C_raw: K*T, temporal components without being deconvolved

% Author: Pengcheng Zhou, Carnegie Mellon University, adapted from Johannes

% options
maxIter = obj.options.maxIter;
if ~exist('smin', 'var') || isempty(smin)
    smin = 3;
end
%% initialization
A = obj.A;
K = size(A, 2);     % number of components
C = obj.C;
C_raw = zeros(size(C));
C_offset = zeros(K, 1);
S = zeros(size(C));
A = full(A);
U = A'*Y;
V = A'*A;
aa = diag(V);   % squares of l2 norm all all components
sn =  zeros(1, K);
kernel = obj.kernel;
kernel_pars = zeros(K, 2);

%% updating
ind_del = false(K, 1);
for miter=1:maxIter
    for k=1:K
        if ind_del
            continue;
        end
        temp = C(k, :) + (U(k, :)-V(k, :)*C)/aa(k);
        % estimate noise
        if miter==1
            sn(k) = get_noise_fft(temp);
        end
        
        % remove baseline
        [temp, C_offset(k)] = remove_baseline(temp, sn(k));
        
        % deconvolution
        if obj.options.deconv_flag
            if miter==1
                [ck, sk, kernel] = deconvCa(temp, kernel, smin, true, false, sn(k));
                kernel_pars(k, :) = kernel.pars;
            else
                kernel.pars = kernel_pars(k, :);
                [ck, sk, kernel] = deconvCa(temp, kernel, smin, false, false, sn(k));
            end
        else
            ck = max(0, temp);
        end
        
        % save convolution kernels and deconvolution results
        C(k, :) = ck;
        
        if sum(ck(2:end))==0
            ind_del(k) = true;
        end
        % save the spike count in the last iteration
        if miter==maxIter
            if obj.options.deconv_flag
                S(k, :) = sk;
            end
            C_raw(k, :) = temp;
        end
    end
end
obj.A = bsxfun(@times, A, sn);
obj.C = bsxfun(@times, C, 1./sn');
obj.C_raw = bsxfun(@times, C_raw, 1./sn');
obj.S = bsxfun(@times, S, 1./sn');
obj.P.kernel_pars = kernel_pars;
obj.delete(ind_del);
