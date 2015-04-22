function [A_or,C_or] = order_ROIs(A,C)

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