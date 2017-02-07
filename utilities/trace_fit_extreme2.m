function max_pr = trace_fit_extreme2(C,YrA,fr,t_int,fac)

if nargin < 5 || isempty(fac); fac = 1; end
if nargin < 4 || isempty(t_int); t_int = 0.25; end
if nargin < 3 || isempty(fr); fr = 30; end

Np = round(t_int*fr);
[K,T] = size(C);
bas = zeros(K,1);
sn = zeros(K,1);
for i = 1:K
    [~,density,xmesh] = kde(YrA(i,:));
    [~,ind] = max(density); 
    bas(i) = xmesh(ind);
    sn(i) = std(YrA(i,YrA(i,:)<bas(i)))/sqrt(1-2/pi);    
end

mu = norminv(1-1/T);
sig = norminv(1-exp(-1)/T) - mu;

z = bsxfun(@times, 1./(fac*sn(:)), bsxfun(@minus, C+YrA, bas+min(C,[],2)));  % normalized z scores
z_pr = exp(-exp(-(z-mu)/sig));

filt_z = filter(ones(1,Np),1,log(z_pr),[],2)/sqrt(Np);
max_pr = exp(max(filt_z,[],2));