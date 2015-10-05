function [A_or,C_or,P_or,srt] = order_ROIs(A,C,P)

% ordering of the found components based on their maximum temporal
% activation and their size (through their l4 norm)

nA = sqrt(sum(A.^2));
nr = length(nA);
A = A/spdiags(nA(:),0,nr,nr);
nA4 = sum(A.^4).^(1/4);
C = spdiags(nA(:),0,nr,nr)*C;
mC = max(C,[],2);
[~,srt] = sort(mC.*nA4','descend');
A_or = A(:,srt);
C_or = C(srt,:);
P_or = P;

if isfield(P,'gn'); P_or.gn=P.gn(srt); end
if isfield(P,'b'); P_or.b=P.b(srt); end
if isfield(P,'c1'); P_or.c1=P.c1(srt); end
if isfield(P,'neuron_sn'); P_or.neuron_sn=P.neuron_sn(srt); end