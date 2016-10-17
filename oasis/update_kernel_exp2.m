function [kernel, s] = update_kernel_exp2(y, s, kernel)
%% estimate convolution kernel with the form of exp(-t/tau_d)-exp(-t/tau_r)

%% inputs:
%       y:  1 X T vector, observed fluorescence trace
%       s:  1 X T vector, spike counts
%       kernel: struct variable with fields {'type', 'fhandle', 'pars',
%       'nMax', 'lb', 'ub'}, it determines the parametric convolution kernel
%           nMax: scalar, length of the convolution kernel
%           lb: 2 X 1 vector, lower bounds of tau_d and tau_r
%           ub: 2 X 1 vector, upper bounds of tau_d and tau_r
%% outputs:
%       kernel: same as input
%       s:  optimal spike train
%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
if ~exist('kernel', 'var') || isempty(kernel)
    kernel = create_kernel('exp2');
elseif ~strcmpi(kernel.type, 'exp2')
    kernel = create_kernel('exp2');
    kernel.nMax = kernel.nMax;
end

bound_pars = kernel.bound_pars;
if bound_pars
    lb = kernel.lb;
    ub = kernel.ub;
else
    lb = kernel.pars/2;
    ub = kernel.pars*2;
end
nMax = kernel.nMax;
fhandle = kernel.fhandle;

t = 1:nMax;
T = length(y);  %number of frames
y = reshape(y, T, 1); % observed fluorescence
s = reshape(s, T, 1); % spike counts

%% create regression matrix
sind = reshape(find(s), 1, []); % indices of all spikes
nspk = length(sind); %  number of vents
if nspk<2 % no events, stop running
    return;
end
temp = bsxfun(@plus, (0:(nMax-1))', sind); % frames that are affected by s
ind = (temp<=T); % remove frames
[ind_tau, ind_s] = find(ind);  %ind_tau corresponds to tau=(t-t'+1); ind_s corresponds to t'; t' is the spike time
yind = temp(ind(:));
temp = sparse(yind, ind_s, ind_tau, T, nspk); % tmtp: t-t'
ind_nonzero = (sum(temp,2)>0); % choose frames affected by s
yv = y(ind_nonzero);
ny = length(yv);    % number of frames used
temp = temp(ind_nonzero, :);
[rsub, csub, tmtp] = find(temp);

%% find the optimal solution by shrinking the searching area
K = 5;
s_all = zeros(K, K, length(sind));
% min_tau1 = lb(1);
% max_tau1 = ub(1);
% min_tau2 = lb(2);
% max_tau2 = ub(2);
f0 = inf;
thresh = 1e-3;
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')
while true
    tau_1 = linspace(lb(1), ub(1), K);
    tau_2 = linspace(lb(2), ub(2), K);
    rss = inf(K);
    
    for m=1:K
        tau_d = tau_1(m);
        for n=1:K
            tau_r = tau_2(n);
            if tau_r>tau_d
                break;
            end
            gt = fhandle([tau_d, tau_r], t);
            H = sparse(rsub, csub, gt(tmtp), ny, nspk);
            sv = (H'*H)\(H'*yv);
            s_all(m, n, :) = sv;
            rss(m, n) = norm(H*sv-yv, 2);
        end
    end
    
    f1 = min(rss(:)); % new residual
    [indr, indc] = find(rss==f1, 1);
    if (f0-f1)/f1 < thresh %improvement is small
        break;
    elseif bound_pars  % shrink the searching area and parameters are bounded
        f0 = f1;
        lb(1) = tau_1(max(1, indr-1));
        ub(1) = tau_1(min(K, indr+1));
        lb(2) = tau_2(max(1, indc-1));
        ub(2) = tau_2(min(K, indc+1));
    else % searching areas are not bounded
        f0 = f1;
        if indr==1
            lb(1) = lb(1)/2;
            ub(1) = tau_1(2);
        elseif indr==K
            lb(1) = tau_1(K-1);
            ub(1) = ub(1)*2;
        else
            lb(1) = tau_1(indr-1);
            ub(1) = tau_1(indr+1);
        end
        
        if indc==1
            lb(2) = lb(2)/2;
            ub(2) = tau_2(2);
        elseif indc==K
            lb(2) = tau_2(K-1);
            ub(2) = ub(2)*2;
        else
            lb(2) = tau_2(indc-1);
            ub(2) = tau_2(indc+1);
        end
        
    end
end
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:singularMatrix')
sv = squeeze(s_all(indr, indc, :));
s(sind) = sv;
kernel.pars = [tau_1(indr), tau_2(indc)];
kernel.fhandle = fhandle;
kernel.nMax = nMax;
end