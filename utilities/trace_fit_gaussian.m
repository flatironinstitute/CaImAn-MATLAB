function fitness = trace_fit_gaussian(C,fr,t_int,fac)

if nargin < 4 || isempty(fac); fac = 3; end
if nargin < 3 || isempty(t_int); t_int = 0.25; end
if nargin < 2 || isempty(fr); fr = 30; end

Np = round(t_int*fr);
[K,T] = size(C);
bas = zeros(K,1);
sn = zeros(K,1);
for i = 1:K
    [~,density,xmesh] = kde(C(i,:));
    [~,ind] = max(density); 
    bas(i) = xmesh(ind);
    sn(i) = std(C(i,C(i,:)<bas(i)))/sqrt(1-2/pi);
end

z = bsxfun(@times, 1./(fac*sn(:)), bsxfun(@minus, C, bas));  % normalized z scores
z_pr = 1 - normcdf(z);

filt_z = filter(ones(1,Np),1,log(z_pr),[],2);
fitness = min(filt_z,[],2);