function [C,f,Y_res,P,S] = update_temporal_components_parallel(Y,A,b,Cin,fin,P,options)

% update temporal components and background given spatial components in
% parallel, by forming of sequence of vertex covers. 
% A variety of different methods can be used and are separated into 2 classes:

% 1-d approaches, where for each component a 1-d trace is computed by removing
% the effect of all the other components and then averaging with the corresponding 
% spatial footprint. Then each trace is denoised. This corresponds to a block-coordinate approach
% 4 different 1-d approaches are included, and any custom method
% can be easily incorporated:
% 'project':                The trace is projected to satisfy the constraints by the (known) calcium indicator dynamics
% 'constrained_foopsi':     The noise constrained deconvolution approach is used. Time constants can be re-estimated (default)
% 'MCEM_foopsi':            Alternating between constrained_foopsi and a MH approach for re-estimating the time constants
% 'MCMC':                   A fully Bayesian method (slowest, but usually most accurate)

% multi-dimensional approaches: (slowest)
% 'noise_constrained': 
% C(j,:) = argmin_{c_j} sum(G*c_j), 
%           subject to:   G*c_j >= 0
%                         ||Y(i,:) - A*C - b*f|| <= sn(i)*sqrt(T)

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
% Y_res:    residual signal (d X T matrix). Y_res = Y - A*C - b*f
% P:        struct for neuron parameters
% S:        deconvolved activity

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

[d,T] = size(Y);
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

if isfield(P,'interp'); Y_interp = P.interp; else Y_interp = sparse(d,T); end        % missing data
if isfield(P,'unsaturatedPix'); unsaturatedPix = P.unsaturatedPix; else unsaturatedPix = 1:d; end   % saturated pixels

mis_data = find(Y_interp);              % interpolate any missing data before deconvolution
Y(mis_data) = Y_interp(mis_data);

if (strcmpi(method,'noise_constrained') || strcmpi(method,'project')) && ~isfield(P,'g')
    options.flag_g = 1;
    if ~isfield(P,'p') || isempty(P.p); P.p = 2; end; 
    p = P.p;
    P = arpfit(Yr,p,options,P.sn);
    if ~iscell(P.g)
        G = make_G_matrix(T,P.g);
    end
else
    G = speye(T);
end

ff = find(sum(A)==0);
if ~isempty(ff)
    A(:,ff) = [];
    if exist('Cin','var')
        if ~isempty(Cin)
            Cin(ff,:) = [];
        end
    end
end

if isempty(Cin) || nargin < 4
    Cin = max((A'*A)\(A'*Y),0);
    ITER = max(ITER,3);
end

if  isempty(b) || isempty(fin) || nargin < 5
    if isempty(b) || nargin < 3
        [b,fin] = nnmf(max(Y - A*Cin,0),1);
    else
        fin = max(b'*Y/norm(b)^2,0);
    end
end

saturatedPix = setdiff(1:d,unsaturatedPix);     % remove any saturated pixels
Ysat = Y(saturatedPix,:);
Asat = A(saturatedPix,:);
bsat = b(saturatedPix,:);
Y = Y(unsaturatedPix,:);
A = A(unsaturatedPix,:);
b = b(unsaturatedPix,:);
d = length(unsaturatedPix);

K = size(A,2);
A = [A,b];
S = zeros(size(Cin));
Cin = [Cin;fin];
C = Cin;

if strcmpi(method,'noise_constrained')
    Y_res = Y - A*Cin;
    mc = min(d,15);  % number of constraints to be considered
    LD = 10*ones(mc,K);
else
    nA = sum(A.^2);
    YrA = Y'*A - Cin'*(A'*A);
    if strcmpi(method,'constrained_foopsi') || strcmpi(method,'MCEM_foopsi')
        P.gn = cell(K,1);
        P.b = cell(K,1);
        P.c1 = cell(K,1);           
        P.neuron_sn = cell(K,1);
    end
    if strcmpi(method,'MCMC')        
        params.B = 300;
        params.Nsamples = 400;
        params.p = P.p;
    else
        params = [];
    end
end
p = P.p;
for iter = 1:ITER
    [O,lo] = update_order(A(:,1:K));
    for jo = 1:length(O)
        Ytemp = YrA(:,O{jo}(:)) + (diag(nA(O{jo}))*Cin(O{jo},:))';
        Ctemp = zeros(length(O{jo}),T);
        Stemp = zeros(length(O{jo}),T);
        btemp = zeros(length(O{jo}),1);
        sntemp = btemp;
        c1temp = btemp;
        gtemp = zeros(length(O{jo}),P.p);
        nT = nA(O{jo});
        parfor jj = 1:length(O{jo})
            if p == 0   % p = 0 (no dynamics assumed)
                cc = max(Ytemp(:,jj)/nT(jj),0);
                Ctemp(jj,:) = full(cc');
                Stemp(jj,:) = C(jj,:);
            else
                switch method
                    case 'project'
                        maxy = max(Ytemp(:,jj)/nT(jj));
                        cc = plain_foopsi(Ytemp(:,jj)/nT(jj)/maxy,G);
                        Ctemp(jj,:) = full(cc')*maxy;
                        Stemp(jj,:) = Ctemp(jj,:)*G';
                    case 'constrained_foopsi'
                        %if restimate_g
                        [cc,cb,c1,gn,sn,spk] = constrained_foopsi(Ytemp(:,jj)/nT(jj),[],[],[],[],options);
                        %else
                        %    [cc,cb,c1,gn,sn,spk] = constrained_foopsi(Ytemp(:,jj)/nA(jj),[],[],P.g,[],options);
                        %end
                        gd = max(roots([1,-gn']));  % decay time constant for initial concentration
                        gd_vec = gd.^((0:T-1));
                        Ctemp(jj,:) = full(cc(:)' + cb + c1*gd_vec);
                        Stemp(jj,:) = spk(:)';
                        Ytemp(:,jj) = Ytemp(:,jj) - nT(jj)*Ctemp(jj,:)';
                        btemp(jj) = cb;
                        c1temp(jj) = c1;
                        sntemp(jj) = sn;   
                        gtemp(jj,:) = gn(:)';
                    case 'MCMC'
                        SAMPLES = cont_ca_sampler(Ytemp(:,jj)/nT(jj),params);
                        Ctemp(jj,:) = make_mean_sample(SAMPLES,Ytemp(:,jj)/nT(jj));
                        Stemp(jj,:) = mean(samples_cell2mat(SAMPLES.ss,T));
                        btemp(jj) = mean(SAMPLES.Cb);
                        c1temp(jj) = mean(SAMPLES.Cin);
                        sntemp(jj) = sqrt(mean(SAMPLES.sn2));
                        gtemp(jj,:) = mean(exp(-1./SAMPLES.g))';
                end
            end
        end
        if p > 0
            if strcmpi(method,'constrained_foopsi') || strcmpi(method,'MCMC');
                P.b(O{jo}) = num2cell(btemp);
                P.c1(O{jo}) = num2cell(c1temp);
                P.neuron_sn(O{jo}) = num2cell(sntemp);
                for jj = 1:length(O{jo})
                    P.gn(O{jo}(jj)) = {gtemp(jj,abs(gtemp(jj,:))>0)'};
                end
                YrA(:,O{jo}(:)) = Ytemp;
                C(O{jo}(:),:) = Ctemp;
                S(O{jo}(:),:) = Stemp;
            end
        end
        fprintf('%i out of %i components updated \n',sum(lo(1:jo)),K);
    end
    ii = K + 1;
    YrA(:,ii) = YrA(:,ii) + nA(ii)*Cin(ii,:)';
    cc = max(YrA(:,ii)/nA(ii),0);
    C(ii,:) = full(cc');
    YrA(:,ii) = YrA(:,ii) - nA(ii)*C(ii,:)';
%     if mod(jj,10) == 0
%         fprintf('%i out of total %i temporal components updated \n',jj,K);
%     end
    %disp(norm(Fin(1:nr,:) - F,'fro')/norm(F,'fro'));
    if norm(Cin - C,'fro')/norm(C,'fro') <= 1e-3
        % stop if the overall temporal component does not change by much
        break;
    else
        Cin = C;
    end
end

if ~strcmpi(method,'noise_constrained')
    Y_res = Y - A*C;
end

f = C(K+1:end,:);
C = C(1:K,:);
Y_res(unsaturatedPix,:) = Y_res;
Y_res(saturatedPix,:) = Ysat - Asat*C - bsat*f;
