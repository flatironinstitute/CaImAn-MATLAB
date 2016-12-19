function [C,f,P,S,YrA] = update_temporal_components_fast(Y,A,b,Cin,fin,P,options)

% update temporal components and background given spatial components
% For each component a 1-d trace is computed by removing
% the effect of all the other components and then averaging with the corresponding 
% spatial footprint. Then each trace is denoised. This corresponds to a block-coordinate approach
% 4 different 1-d approaches are included, and any custom method
% can be easily incorporated:
% 'constrained_foopsi':     The noise constrained deconvolution approach is used. Time constants can be re-estimated (default)
% 'MCMC':                   A fully Bayesian method (slowest, but usually most accurate)

% INPUTS:
% Y:        raw data ( d X T matrix)
% A:        spatial footprints  (d x nr matrix)
% b:        spatial background  (d x 1 vector)
% Cin:      current estimate of temporal components (nr X T matrix)
% fin:      current estimate of temporal background (1 x T vector)
% P:        struct for neuron parameters
% options:  struct for algorithm parameters

% LD:       Lagrange multipliers (needed only for 'noise_constrained' method).
% 
% OUTPUTS:
% C:        temporal components (nr X T matrix)
% f:        temporal background (1 x T vector)
% P:        struct for neuron parameters
% S:        deconvolved activity

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2016

memmaped = isobject(Y);
if memmaped
    sizY = size(Y,'Y');
    d = prod(sizY(1:end-1));
    T = sizY(end);
else
    [d,T] = size(Y);
end
if isempty(P) || nargin < 6
    active_pixels = find(sum(A,2));                                 % pixels where the greedy method found activity
    unsaturated_pixels = find_unsaturatedPixels(Y);                 % pixels that do not exhibit saturation
    options.pixels = intersect(active_pixels,unsaturated_pixels);   % base estimates only on unsaturated, active pixels                
end

defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = []; end
if ~isfield(options,'deconv_method') || isempty(options.deconv_method); method = defoptions.deconv_method; else method = options.deconv_method; end  % choose method
if ~isfield(options,'restimate_g') || isempty(options.restimate_g); restimate_g = defoptions.restimate_g; else restimate_g = options.restimate_g; end % re-estimate time constant (only with constrained foopsi)
if ~isfield(options,'temporal_iter') || isempty(options.temporal_iter); ITER = defoptions.temporal_iter; else ITER = options.temporal_iter; end           % number of block-coordinate descent iterations
if ~isfield(options,'bas_nonneg'); options.bas_nonneg = defoptions.bas_nonneg; end
if ~isfield(options,'fudge_factor'); options.fudge_factor = defoptions.fudge_factor; end

if isfield(P,'unsaturatedPix'); unsaturatedPix = P.unsaturatedPix; else unsaturatedPix = 1:d; end   % saturated pixels

ff = find(sum(A)==0);
if ~isempty(ff)
    A(:,ff) = [];
    if exist('Cin','var')
        if ~isempty(Cin)
            Cin(ff,:) = [];
        end
    end
end

% estimate temporal (and spatial) background if they are not present
if isempty(fin) || nargin < 5           % temporal background missing
    if isempty(b) || nargin < 3
        fin = mm_fun(ones(d,1),Y);
        fin = fin/norm(fin);
        b = max(mm_fun(fin,Y),0);
        options.nb = 1;
    else
        fin = max(b(bk_pix,:)'*Y(bk_pix,:),0)/(b(bk_pix,:)'*b(bk_pix,:));
    end
end

% construct product A'*Y

K = size(A,2);
nb = size(b,2);

AY = mm_fun([A,double(b)],Y);
bY = AY(K+1:end,:);
AY = AY(1:K,:);

if isempty(Cin) || nargin < 4    % estimate temporal components if missing    
    Cin = max((A'*A)\double(AY - (A'*b)*fin),0);  
end

A = [A,double(b)];
S = zeros(size(Cin));
Cin = [Cin;fin];
C = Cin;

AA = (A'*A);
AY = [AY;bY];

if strcmpi(method,'constrained_foopsi')
    P.gn = cell(K,1);
    P.b = num2cell(zeros(K,1));
    P.c1 = num2cell(zeros(K,1));           
    P.neuron_sn = num2cell(zeros(K,1));
end
if strcmpi(method,'MCMC')        
    params.B = 300;
    params.Nsamples = 400;
    params.p = P.p;
    params.bas_nonneg = options.bas_nonneg;
else
    params = [];
end

p = P.p;
options.p = P.p;
C = double(C);

C = HALS_temporal(AY, A, C, 100, [], options.bas_nonneg, true);
if p > 0
    YrA = bsxfun(@times, 1./sum(A.^2)',AY - AA*C);
    if options.temporal_parallel
        C = mat2cell(C,ones(size(C,1),1),T);
        YrA = mat2cell(YrA,ones(size(C,1),1),T);
        S = cell(K,1);   
        b = cell(K,1);
        c1 = cell(K,1);
        gn = cell(K,1);
        neuron_sn = cell(K,1);
        %disp([K,nb,length(C)])
        use_OASIS = false;
        parfor ii = 1:K+nb
            Ytemp = C{ii} + YrA{ii};
            if ii <= K
                %fprintf(num2str(ii))
                if strcmpi(method,'MCMC')
                    samples_mcmc = deal(struct('lam_', [], 'spiketimes_', [], 'A_', [], 'b_', [], 'C_in', [], 'sg', [], 'g', []));
                end                 
                switch method
                    case 'constrained_foopsi'
                        [cc,b_temp,c1_temp,gn_temp,neuron_sn_temp,spk] = constrained_foopsi(Ytemp,[],[],[],[],options);
                            %P.gn{ii} = gn;
                        gd = max(roots([1,-gn_temp']));  % decay time constant for initial concentration
                        gd_vec = gd.^((0:T-1));
                        C{ii} = full(cc(:)' + b_temp + c1_temp*gd_vec);
                        S{ii} = spk(:)';
                        b{ii} = b_temp;
                        c1{ii} = c1_temp;
                        neuron_sn{ii} = neuron_sn_temp;
                        gn{ii} = gn_temp;
                        if use_OASIS
                            [cc,spk,kernel] = deconvCa(Ytemp,[],[],true,false,[],5);                        
                            S{ii} = spk(:)';
                            b_temp = median(Ytemp-cc(:)')
                            C{ii} = cc(:)' + b_temp;
                            b{ii} = b_temp;
                            c1{ii} = Ytemp(1)-cc(1);
                            neuron_sn{ii} = std(Ytemp-cc(:)');
                            gn{ii} = kernel.pars;                        
                        else
                            [cc,b_temp,c1_temp,gn_temp,neuron_sn_temp,spk] = constrained_foopsi(Ytemp,[],[],[],[],options);
                            %P.gn{ii} = gn;
                            gd = max(roots([1,-gn_temp']));  % decay time constant for initial concentration
                            gd_vec = gd.^((0:T-1));
                            C{ii} = full(cc(:)' + b_temp + c1_temp*gd_vec);
                            S{ii} = spk(:)';
                            b{ii} = b_temp;
                            c1{ii} = c1_temp;
                            neuron_sn{ii} = neuron_sn_temp;
                            gn{ii} = gn_temp;
                        end
                    case 'MCMC'
                        SAMPLES = cont_ca_sampler(Ytemp,params);
                        ctemp = make_mean_sample(SAMPLES,Ytemp);
                        C{ii} = ctemp(:)';
                        S{ii} = mean(samples_cell2mat(SAMPLES.ss,T));
                        b{ii} = mean(SAMPLES.Cb);
                        c1{ii} = mean(SAMPLES.Cin);
                        neuron_sn{ii} = sqrt(mean(SAMPLES.sn2));
                        gr = mean(exp(-1./SAMPLES.g));
                        gp = poly(gr);
                        gn{ii} = -gp(2:end);
                        samples_mcmc(ii) = SAMPLES; % FN added, a useful parameter to have.
                end                
            else
                C{ii} = max(Ytemp(:),0)';
            end
        end
        C = cell2mat(C);
        S = cell2mat(S);
        %YrA = cell2mat(YrA);
        P.b = b;
        P.c1 = c1;
        P.neuron_sn = neuron_sn;
        P.gn = gn;
    else
        for ii = 1:K+nb
            Ytemp = C(ii,:) + YrA(ii,:);
            if ii <= K
                switch method
                     case 'constrained_foopsi'
%                         if restimate_g
%                             [cc,cb,c1,gn,sn,spk] = constrained_foopsi(Ytemp,[],[],[],[],options);
%                             P.gn{ii} = gn;
%                         else
%                             [cc,cb,c1,gn,sn,spk] = constrained_foopsi(Ytemp,[],[],P.g,[],options);
%                         end
%                         gd = max(roots([1,-gn']));  % decay time constant for initial concentration
%                         gd_vec = gd.^((0:T-1));
%                         C(ii,:) = full(cc(:)' + cb + c1*gd_vec);
%                         S(ii,:) = spk(:)';
%                         P.b{ii} = cb;
%                         P.c1{ii} = c1;           
%                         P.neuron_sn{ii} = sn;
                        [cc,spk,kernel] = deconvCa(Ytemp,[],[],true,false,[],5);                        
                        S(ii,:) = spk(:)';
                        P.b{ii} = median(Ytemp-cc(:)');
                        C(ii,:) = cc(:)'+P.b{ii};
                        P.c1{ii} = Ytemp(1)-cc(1);
                        P.neuron_sn{ii} = std(Ytemp-cc(:)');
                        P.gn{ii} = kernel.pars;                            
                    case 'MCMC'
                        SAMPLES = cont_ca_sampler(Ytemp,params);
                        ctemp = make_mean_sample(SAMPLES,Ytemp);
                        C(ii,:) = ctemp';
                        S(ii,:) = mean(samples_cell2mat(SAMPLES.ss,T));
                        P.b{ii} = mean(SAMPLES.Cb);
                        P.c1{ii} = mean(SAMPLES.Cin);
                        P.neuron_sn{ii} = sqrt(mean(SAMPLES.sn2));
                        gr = mean(exp(-1./SAMPLES.g));
                        gp = poly(gr);
                        P.gn{ii} = -gp(2:end);
                        P.samples_mcmc(ii) = SAMPLES; % FN added, a useful parameter to have.
                end            
            else
                C(ii,:) = max(Ytemp(:),0)';
            end
        end
    end
else
    warning('No deconvolution is performed. \n');
end

YrA = bsxfun(@times, 1./sum(A.^2)',AY - AA*C);
f = C(K+1:end,:);
C = C(1:K,:);
YrA = YrA(1:K,:);