function [C, P, S] = deconv_temporal(C,P,options)
% deconvolve all temporal components independently using FOOPSI
% The noise constrained deconvolution approach is used. Time constants can be re-estimated (default)

% The update can happen either in parallel (default) or serial by tuning options.temporal_parallel.
% In the case of parallel implementation the methods 'MCEM_foopsi' and 'noise_constrained' are not supported

% INPUTS:
% C:      current estimate of temporal components (nr X T matrix)
% P:        struct for neuron parameters
% options:  struct for algorithm parameters

% OUTPUTS:
% C:        temporal components (nr X T matrix)
% P:        struct for neuron parameters
% S:        deconvolved activity

% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
% it's modified from udpate_temporal_components.m written by Eftychios

defoptions = CNMFSetParms;
if nargin<2;    P.p = 2; end
if P.p==0; return; end

if ~isfield(options,'restimate_g') || isempty(options.restimate_g); restimate_g = defoptions.restimate_g; else restimate_g = options.restimate_g; end % re-estimate time constant (only with constrained foopsi)
if ~isfield(options,'bas_nonneg'); options.bas_nonneg = defoptions.bas_nonneg; end
if ~isfield(options,'fudge_factor'); options.fudge_factor = defoptions.fudge_factor; end
if ~isfield(options,'temporal_parallel'); options.temporal_parallel = defoptions.temporal_parallel; end

options.p = P.p;

[K, T] = size(C);

if ~isfield(P, 'gn'); options.restimate=false; end
if options.restimate_g; P.gn = cell(K,1); end
P.b = cell(K,1);
P.c1 = cell(K,1);
P.neuron_sn = cell(K,1);

S = zeros(K, T);
Cb = zeros(K,1);
Sn = zeros(K,1);


fprintf('Progress:\n');
fprintf(['\n' repmat('.',1,K) '\n\n']);
if options.temporal_parallel
    tmp_gn = zeros(K, P.p);
    if ~restimate_g
        for jj=1:K
            tmp_gn(jj, :) =  P.gn{jj};
        end
    end
    parfor jj = 1:K
        if restimate_g
            [cc,cb,c1,gn,sn,spk] = constrained_foopsi(C(jj,:),[],[],[],[],options);
        else
            [cc,cb,c1,gn,sn,spk] = constrained_foopsi(C(jj,:),[],[],[],tmp_gn(jj,:),options);
        end
        tmp_gn(jj, :) = gn;
        gd = max(roots([1,-gn']));  % decay time constant for initial concentration
        gd_vec = gd.^((0:T-1));
        C(jj,:) = full(cc(:)'+ c1*gd_vec)+cb;
        S(jj,:) = spk(:)';
        Sn(jj) = sn;
        
        fprintf('\b|\n');
    end
    for jj=1:K
        P.gn{jj} = tmp_gn(jj, :);
    end
else
    for jj = 1:K
        if restimate_g
            [cc,cb,c1,gn,sn,spk] = constrained_foopsi(C(jj,:),[],[],[],[],options);
        else
            [cc,cb,c1,gn,sn,spk] = constrained_foopsi(C(jj,:),[],[],[],P.gn{jj},options);
        end
        P.gn{jj} = gn;
        gd = max(roots([1,-gn']));  % decay time constant for initial concentration
        gd_vec = gd.^((0:T-1));
        C(jj,:) = full(cc(:)'+ c1*gd_vec)+cb;
        S(jj,:) = spk(:)';
        Sn(jj) = sn;
        
        fprintf('\b|\n');
    end
end
P.neuron_sn = Sn;