function [A_or,C_or,S_or,P_or,srt] = order_ROIs(A,C,S,P, srt)

% ordering of the found components based on their maximum temporal
% activation and their size (through their l_inf norm)
% you can also pre-specify the ordering sequence

nA = sqrt(sum(A.^2));
nr = length(nA);
A = A/spdiags(nA(:),0,nr,nr);
C = spdiags(nA(:),0,nr,nr)*C;
mA = sum(A.^4).^(1/4);
%sA = sum(A);
mC = max(C,[],2);
if ~exist('srt', 'var')||isempty(srt)
    [~,srt] = sort(mC.*mA','descend');
end
A_or = A(:,srt);
C_or = C(srt,:);

if nargin < 4
    P_or = [];
else
    P_or = P;
    if isfield(P,'gn'); P_or.gn=P.gn(srt); end
    if isfield(P,'b'); P_or.b=P.b(srt); end
    if isfield(P,'c1'); P_or.c1=P.c1(srt); end
    if isfield(P,'neuron_sn'); P_or.neuron_sn=P.neuron_sn(srt); end
end

if nargin < 3 || isempty(S)
    S_or = [];
else
    S = spdiags(nA(:),0,nr,nr)*S;
    S_or = S(srt,:);
end