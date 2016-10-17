function SAMPLES = cont_ca_sampler(Y,params)

% MCMC continuous time sampler of spiketimes given a fluorescence trace Y

% Inputs:
% Y                     data
% params                parameters structure
% params.g              discrete time constant(s) (estimated if not provided)
% params.sn             initializer for noise (estimated if not provided)
% params.b              initializer for baseline (estimated if not provided)
% params.c1             initializer for initial concentration (estimated if not provided)
% params.Nsamples       number of samples after burn in (default 500)
% params.B              number of burn in samples (default 200)
% params.marg           flag for marginalized sampler for baseline and initial concentration (default 0)
% params.upd_gam        flag for updating gamma (default 1)
% params.gam_step       number of samples after which gamma is updated (default 50)
% params.std_move       standard deviation of shifting kernel (default 3*Dt)
% params.add_move       number of add moves per iteration (default T/100)
% params.init           initial sample 
% params.f              imaging rate (default 1)
% params.p              order of AR model (p == 1 or p == 2, default 1)
% params.defg           default discrete time constants in case constrained_foopsi cannot find stable estimates
% params.TauStd         standard deviation for time constants in continuous time (default [0.2,2])
% params.A_lb           lower bound for spike amplitude (default 0.1*range(Y))
% params.b_lb           lower bound for baseline (default 0.01 quantile of Y)
% params.c1_lb          lower bound for initial value (default 0)

% output struct SAMPLES
% ss                    Nsamples x 1 cells with spike times for each sample
% ns                    Nsamples x 1 vector with number of spikes
% Am                    Nsamples x 1 vector with samples for spike amplitude
% ld                    Nsamples x 1 vector with samples for firing rate

% If marginalized sampler is used (params.marg = 1)
% Cb                    posterior mean and sd for baseline
% Cin                   posterior mean and sd for initial condition
% else
% Cb                    Nsamples x 1 vector with samples for baseline
% Cin                   Nsamples x 1 vector with samples for initial concentration
% sn                    Nsamples x 1 vector with samples for noise variance

% If gamma is updated (params.upd_gam = 1)
% g                     Nsamples x p vector with the gamma updates

% Author: Eftychios A. Pnevmatikakis and Josh Merel

% Reference: Pnevmatikakis et al., Bayesian spike inference from calcium
%                imaging data, Asilomar 2013.

Y = Y(:);
T = length(Y);
isanY = ~isnan(Y);                          % Deal with possible missing entries
E = speye(T);
E = E(isanY,:);

% define default parameters
defparams.g = [];                           % initializer for time constants
defparams.sn = [];                          % initializer for noise std
defparams.b = [];                           % initializer for baseline concentration
defparams.c1 = [];                          % initializer for initial concentration
defparams.c = [];                           % initializer for calcium concentration
defparams.sp = [];                          % initializer for spiking signal
defparams.bas_nonneg = 0;                   % allow negative baseline during initialization
defparams.Nsamples = 400;                   % number of samples after burn in period
defparams.B = 200;                          % length of burn in period
defparams.marg = 0;                         % flag to marginalize out baseline and initial concentration
defparams.upd_gam = 1;                      % flag for updating time constants
defparams.gam_step = 1;                     % flag for how often to update time constants
defparams.A_lb = 0.1*range(Y);              % lower bound for spike amplitude
defparams.b_lb = quantile(Y,0.01);          % lower bound for baseline
defparams.c1_lb = 0;                        % lower bound for initial concentration

defparams.std_move = 3;                     % standard deviation of spike move kernel
defparams.add_move = ceil(T/100);           % number of add moves
defparams.init = [];                        % sampler initializer
defparams.f = 1;                            % imaging rate (irrelevant)
defparams.p = 1;                            % order of AR process (use p = 1 or p = 2)
defparams.defg = [0.6,0.95];                % default time constant roots
defparams.TauStd = [.2,2];                  % Standard deviation for time constant proposal
defparams.prec = 1e-2;                      % Precision parameter when adding new spikes
defparams.con_lam = true;                   % Flag for constant firing across time
defparams.print_flag = 0;

if nargin < 2
    params = defparams;
else
    if ~isfield(params,'g'); params.g = defparams.g; end
    if ~isfield(params,'sn'); params.sn = defparams.sn; end
    if ~isfield(params,'b'); params.b = defparams.b; end
    if ~isfield(params,'c1'); params.c1 = defparams.c1; end
    if ~isfield(params,'c'); params.c = defparams.c; end
    if ~isfield(params,'sp'); params.sp = defparams.sp; end
    if ~isfield(params,'bas_nonneg'); params.bas_nonneg = defparams.bas_nonneg; end
    if ~isfield(params,'Nsamples'); params.Nsamples = defparams.Nsamples; end
    if ~isfield(params,'B'); params.B = defparams.B; end
    if ~isfield(params,'marg'); params.marg = defparams.marg; end
    if ~isfield(params,'upd_gam'); params.upd_gam = defparams.upd_gam; end
    if ~isfield(params,'gam_step'); params.gam_step = defparams.gam_step; end
    if ~isfield(params,'std_move'); params.std_move = defparams.std_move; end
    if ~isfield(params,'add_move'); params.add_move = defparams.add_move; end
    if ~isfield(params,'init'); params.init = defparams.init; end
    if ~isfield(params,'f'); params.f = defparams.f; end
    if ~isfield(params,'p'); params.p = defparams.p; end
    if ~isfield(params,'defg'); params.defg = defparams.defg; end
    if ~isfield(params,'TauStd'); params.TauStd = defparams.TauStd; end
    if ~isfield(params,'A_lb'); params.A_lb = defparams.A_lb; end
    if ~isfield(params,'b_lb'); params.b_lb = defparams.b_lb; end
    if ~isfield(params,'c1_lb'); params.c1_lb = defparams.c1_lb; end
    if ~isfield(params,'prec'); params.prec = defparams.prec; end
    if ~isfield(params,'con_lam'); params.con_lam = defparams.con_lam; end
    if ~isfield(params,'print_flag'); params.print_flag = defparams.print_flag; end    
end
       
Dt = 1;                                     % length of time bin

marg_flag = params.marg;
gam_flag = params.upd_gam;
gam_step = params.gam_step;
std_move = params.std_move;
add_move = params.add_move;

if isempty(params.init)
   params.init = get_initial_sample(Y,params);   
end 
SAM = params.init;

if isempty(params.g)
    p = params.p;
else
    p = length(params.g);                       % order of autoregressive process
end

g = SAM.g(:)';   % check initial time constants, if not reasonable set to default values 

if g == 0 
    gr = params.defg;
    pl = poly(gr);
    g = -pl(2:end);
    p = 2;
end
gr = sort(roots([1,-g(:)']));
if p == 1; gr = [0,max(gr)]; end
if any(gr<0) || any(~isreal(gr)) || length(gr)>2 || max(gr)>0.998
    gr = params.defg;
end
tau = -Dt./log(gr);
tau1_std = max(tau(1)/100,params.TauStd(1));
tau2_std = min(tau(2)/5,params.TauStd(2)); 
ge = max(gr).^(0:T-1)';
if p == 1
    G1 = sparse(1:T,1:T,Inf*ones(T,1));
elseif p == 2
    G1 = spdiags(ones(T,1)*[-min(gr),1],[-1:0],T,T);
else
    error('This order of the AR process is currently not supported');
end
G2 = spdiags(ones(T,1)*[-max(gr),1],[-1:0],T,T);

sg = SAM.sg;

spiketimes_ = SAM.spiketimes_;
lam_ = SAM.lam_;
A_ = SAM.A_*diff(gr);
b_ = max(SAM.b_,prctile(Y,8));
C_in = max(min(SAM.C_in,Y(1)-b_),0);
    
s_1 = zeros(T,1);
s_2 = zeros(T,1);
s_1(ceil(spiketimes_/Dt)) = exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau(1));
s_2(ceil(spiketimes_/Dt)) = exp((spiketimes_ - Dt*ceil(spiketimes_/Dt))/tau(2));    

prec = params.prec;

ef_d = exp(-(0:T)/tau(2));   % construct transient exponentials
if p == 1
    h_max = 1;              % max value of transient    
    ef_h = [0,0];
    e_support = find(ef_d<prec*h_max,1,'first');
    if isempty(e_support);
        e_support = T;
    end
    e_support = min(e_support,T);
else
    t_max = (tau(1)*tau(2))/(tau(2)-tau(1))*log(tau(2)/tau(1)); %time of maximum
    h_max = exp(-t_max/tau(2)) - exp(-t_max/tau(1));            % max value of transient    
    ef_h = -exp(-(0:T)/tau(1));
    e_support = find(ef_d-ef_h<prec*h_max,1,'first');
    if isempty(e_support);
        e_support = T;
    end
    e_support = min(e_support,T);
end
ef_h = ef_h(1:min(e_support,length(ef_h)))/diff(gr);
ef_d = ef_d(1:e_support)/diff(gr);
ef = [{ef_h ef_d};{cumsum(ef_h.^2) cumsum(ef_d.^2)}];


B = params.B;
N = params.Nsamples + B;
if p == 1; G1sp = zeros(T,1); else G1sp = G1\s_1(:); end
Gs = (-G1sp(:)+G2\s_2(:))/diff(gr);

ss = cell(N,1); 
lam = zeros(N,1);
Am = zeros(N,1);
ns = zeros(N,1);
Gam = zeros(N,2);
if ~marg_flag
    Cb = zeros(N,1);
    Cin = zeros(N,1);
    SG = zeros(N,1);
end

Sp = .1*range(Y)*eye(3);                          % prior covariance for [A,Cb,Cin]
Ld = inv(Sp);
lb = [params.A_lb/h_max*diff(gr),params.b_lb,params.c1_lb]';      % lower bound for [A,Cb,Cin]
A_ = max(A_,1.1*lb(1));

mu = [A_;b_;C_in];                % prior mean 
Ns = 15;                          % Number of HMC samples

Ym = Y - ones(T,1)*mu(2) - ge*mu(3);

mub = zeros(2,1);
Sigb = zeros(2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra tau-related params
%%%%%%%%%%%%%%%%%%%%%%%%%%%

tauMoves = [0 0];
tau_min = 0;
tau_max = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:N 
    if gam_flag
        Gam(i,:) = tau;
    end
    sg_ = sg;
    rate = @(t) lambda_rate(t,lam_);
    [spiketimes, ~]  = get_next_spikes(spiketimes_(:)',A_*Gs',Ym',ef,tau,sg_^2, rate, std_move, add_move, Dt, A_, params.con_lam);
    spiketimes(spiketimes<0) = -spiketimes(spiketimes<0);
    spiketimes(spiketimes>T*Dt) = 2*T*Dt - spiketimes(spiketimes>T*Dt); 
    spiketimes_ = spiketimes;    
    trunc_spikes = ceil(spiketimes/Dt);
    trunc_spikes(trunc_spikes == 0) = 1;
    s_1 = zeros(T,1);
    s_2 = zeros(T,1);
    s_1(trunc_spikes) = exp((spiketimes_ - Dt*trunc_spikes)/tau(1));
    s_2(trunc_spikes) = exp((spiketimes_ - Dt*trunc_spikes)/tau(2));   
    if p == 1; G1sp = zeros(T,1); else G1sp = G1\s_1(:); end
    Gs = (-G1sp+G2\s_2(:))/diff(gr);
    
    %s_ = s_2 - s_1 - min(gr)*[0;s_2(1:end-1)] + max(gr)*[0;s_1(1:end-1)];
    %Gs = filter(1,[1,-sum(gr),prod(gr)]',s_/diff(gr));
    
    ss{i} = spiketimes;
    nsp = length(spiketimes);
    ns(i) = nsp;
    lam(i) = nsp/(T*Dt);
    lam_ = lam(i);
    
    AM = [Gs,ones(T,1),ge];
    EAM = E*AM;
    L = inv(Ld + EAM'*EAM/sg^2);
    mu_post = (Ld + EAM'*EAM/sg^2)\(EAM'*Y(isanY)/sg^2 + Sp\mu);
    if ~marg_flag
        x_in = [A_;b_;C_in];
        if any(x_in <= lb)
            x_in = max(x_in,(1+0.1*sign(lb)).*lb) + 1e-5;
        end
        if all(isnan(L(:))) % FN added to avoid error in R = chol(L) in HMC_exact2 due to L not being positive definite. It happens when isnan(det(Ld + AM'*AM/sg^2)), ie when Ld + AM'*AM/sg^2 is singular (not invertible).
            Am(i) = NaN;
            Cb(i) = NaN;
            Cin(i) = NaN';
        else
            [temp,~] = HMC_exact2(eye(3), -lb, L, mu_post, 1, Ns, x_in);
            Am(i) = temp(1,Ns);
            Cb(i) = temp(2,Ns);
            Cin(i) = temp(3,Ns)';
        end
        A_ = Am(i);
        b_ = Cb(i);
        C_in = Cin(i);

        Ym   = Y - b_ - ge*C_in;
        res   = Ym - A_*Gs;
        sg   = 1./sqrt(gamrnd(1+length(isanY)/2,1/(0.1 + sum((res(isanY).^2)/2))));
        SG(i) = sg;
    else
        repeat = 1;
        cnt = 0;
        while repeat
            A_ = mu_post(1) + sqrt(L(1,1))*randn;
            repeat = (A_<0);
            cnt = cnt + 1;
            if cnt > 1e3
                error('The marginalized sampler cannot find a valid amplitude value. Set params.marg = 0 and try again.')
            end
        end                
        Am(i) = A_;
        if i > B
           mub = mub + mu_post(2:3); %mu_post(1+(1:p));
           Sigb = Sigb + L(2:3,2:3); %L(1+(1:p),1+(1:p));
        end
    end
    if gam_flag
        if mod(i-B,gam_step) == 0  % update time constants
            if p >= 2       % update rise time constant
                logC = -norm(Y(isanY) - EAM*[A_;b_;C_in])^2; 
                tau_ = tau;
                tau_temp = tau_(1)+(tau1_std*randn); 
                while tau_temp >tau(2) || tau_temp<tau_min
                    tau_temp = tau_(1)+(tau1_std*randn);
                end 
                tau_(1) = tau_temp;
                gr_ = exp(Dt*(-1./tau_));
                s_1_ = zeros(T,1);
                s_1_(trunc_spikes) = exp((spiketimes_ - Dt*trunc_spikes)/tau_(1));
                G1_ = spdiags(ones(T,1)*[-min(gr_),1],[-1:0],T,T);
                Gs_ = (-G1_\s_1_(:)+G2\s_2(:))/diff(gr_);

                logC_ = -norm(E*(Y(:)-A_*Gs-b_-C_in*ge))^2;
                %accept or reject
                prior_ratio = 1;
                ratio = exp((logC_-logC)/(2*sg^2))*prior_ratio;
                if rand < ratio %accept
                    tau = tau_;
                    G1 = G1_; %c = c_; 
                    Gs = Gs_; gr = gr_; s_1 = s_1_;
                    tauMoves = tauMoves + [1 1];
                else
                    tauMoves = tauMoves + [0 1];
                end                
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            % next update decay time constant
            %%%%%%%%%%%%%%%%%%%%%%%

            %initial logC
            logC = -norm(E*(Y(:)-A_*Gs-b_-C_in*ge))^2;
            tau_ = tau;
            tau_temp = tau_(2)+(tau2_std*randn);
            while tau_temp>tau_max || tau_temp<tau_(1)
                tau_temp = tau_(2)+(tau2_std*randn);
            end  
            tau_(2) = tau_temp;
            s_2_ = zeros(T,1);
            s_2_(trunc_spikes) = exp((spiketimes_ - Dt*trunc_spikes)/tau_(2));
            gr_ = exp(Dt*(-1./tau_));
            ge_ = max(gr_).^(0:T-1)';
            if p == 1; G1sp = zeros(T,1); else G1sp = G1\s_1(:); end
            G2_ = spdiags(ones(T,1)*[-max(gr_),1],[-1:0],T,T);              
            Gs_ = (-G1sp+G2_\s_2_(:))/diff(gr_);

            logC_ = -norm(E*(Y(:)-A_*Gs_-b_-C_in*ge_))^2;

            %accept or reject
            prior_ratio = 1;
            ratio = exp((1./(2*sg^2)).*(logC_-logC))*prior_ratio;
            if rand<ratio %accept
                tau = tau_;
                %c = c_; 
                ge = ge_; Gs = Gs_; G2 = G2_; gr = gr_; s_2 = s_2_;
                tauMoves = tauMoves + [1 1];
            else
                tauMoves = tauMoves + [0 1];
            end
            
            ef_d = exp(-(0:T)/tau(2));
            if p == 1
                h_max = 1;              % max value of transient    
                ef_h = [0,0];
                e_support = find(ef_d<prec*h_max,1,'first');
                if isempty(e_support);
                    e_support = T;
                end
                e_support = min(e_support,T);
            else
                t_max = (tau(1)*tau(2))/(tau(2)-tau(1))*log(tau(2)/tau(1)); %time of maximum
                h_max = exp(-t_max/tau(2)) - exp(-t_max/tau(1));            % max value of transient    
                ef_h = -exp(-(0:T)/tau(1));
                e_support = find(ef_d-ef_h<prec*h_max,1,'first');
                if isempty(e_support);
                    e_support = T;
                end
                e_support = min(e_support,T);
            end
            ef_h = ef_h(1:min(e_support,length(ef_h)))/diff(gr);
            ef_d = ef_d(1:e_support)/diff(gr);
            ef = [{ef_h ef_d};{cumsum(ef_h.^2) cumsum(ef_d.^2)}];     
        end
    end
    if params.print_flag && mod(i,100)==0
        fprintf('%i out of total %i samples drawn \n', i, N);
    end 
end
if marg_flag
    mub = mub/(N-B);
    Sigb = Sigb/(N-B)^2;
end

if marg_flag
    SAMPLES.Cb = [mub(1),sqrt(Sigb(1,1))];
    SAMPLES.Cin = [mub(2),sqrt(Sigb(2,2))]; %[mub(1+(1:p)),sqrt(diag(Sigb(1+(1:p),1+(1:p))))];
else
    SAMPLES.Cb = Cb(B+1:N);
    SAMPLES.Cin = Cin(B+1:N,:);
    SAMPLES.sn2 = SG(B+1:N).^2;
end
SAMPLES.ns = ns(B+1:N);
SAMPLES.ss = ss(B+1:N);
SAMPLES.ld = lam(B+1:N);
SAMPLES.Am = Am(B+1:N);
if gam_flag
    SAMPLES.g = Gam(B+1:N,:);
end
SAMPLES.params = params.init;