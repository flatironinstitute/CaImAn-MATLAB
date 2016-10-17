function c_m = make_mean_sample(SAMPLES,Y)

% construct mean calcium sample from samples structure SAMPLES

% Inputs:
% SAMPLES:      samples struct output of cont_ca_sampler
% Y:            observed fluorescence trace

% Output:
% c_m:          mean calcium sample

% Author: Eftychios A. Pnevmatikakis, 2016, Simons Foundation.

T = length(Y);
N = length(SAMPLES.ns);
P = SAMPLES.params;
P.f = 1;
g = P.g(:);
Dt = 1/P.f;                                     % length of time bin
if ~isfield(SAMPLES,'g');
    SAMPLES.g = ones(N,1)*g';
end


if length(SAMPLES.Cb) == 2
    marg = 1;       % marginalized sampler
else
    marg = 0;       % full sampler
end

C_rec = zeros(N,T);
for rep = 1:N
    %trunc_spikes = ceil(SAMPLES.ss{rep}/Dt);
    tau = SAMPLES.g(rep,:);
    gr = exp(-1./tau);

    ge = max(gr).^(0:T-1)';
    s_1 =   sparse(ceil(SAMPLES.ss{rep}/Dt),1,exp((SAMPLES.ss{rep} - Dt*ceil(SAMPLES.ss{rep}/Dt))/tau(1)),T,1);  
    s_2 =   sparse(ceil(SAMPLES.ss{rep}/Dt),1,exp((SAMPLES.ss{rep} - Dt*ceil(SAMPLES.ss{rep}/Dt))/tau(2)),T,1);  
    if gr(1) == 0
        G1 = sparse(1:T,1:T,Inf*ones(T,1));
        G1sp = zeros(T,1);
    else
        G1 = spdiags(ones(T,1)*[-min(gr),1],[-1:0],T,T);
        G1sp = G1\s_1(:);
    end    
    G2 = spdiags(ones(T,1)*[-max(gr),1],[-1:0],T,T);
    Gs = (-G1sp+ G2\s_2(:))/diff(gr);
    if marg
        %C_rec(rep,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*filter(1,[1,-SAMPLES.g(rep,:)],full(s_)+[SAMPLES.Cin(:,1)',zeros(1,T-p)]);
        C_rec(rep,:) = SAMPLES.Cb(1) + SAMPLES.Am(rep)*Gs + (ge*SAMPLES.Cin(:,1));
    else
        %C_rec(rep,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*filter(1,[1,-SAMPLES.g(rep,:)],full(s_)+[SAMPLES.Cin(rep,:),zeros(1,T-p)]);
        C_rec(rep,:) = SAMPLES.Cb(rep) + SAMPLES.Am(rep)*Gs + (ge*SAMPLES.Cin(rep,:)');
    end
end

c_m = mean(C_rec);