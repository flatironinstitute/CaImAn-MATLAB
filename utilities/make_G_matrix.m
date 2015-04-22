function G = make_G_matrix(T,g,varargin)

if length(g) == 1 && g < 0
    g=0;
end

G = spdiags(ones(T,1)*[-flipud(g(:))',1],-length(g):0,T,T);
if nargin == 3
    sl = [0;cumsum(varargin{1}(:))];
    for i = 1:length(sl)-1
        G(sl(i)+1,sl(i+1))=0;
    end
end