function [C_df,Df] = extract_DF_F(Y,A,C,i)

% extract DF/F signals after performing NMF
% inputs: Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%         A matrix of spatial components (d x K matrix, K # of components)
%         C matrix of temporal components (K x T matrix)
%         index of component that represent the background (optional, if not
%         given it's estimated)

% outputs:  C_df temporal components in the DF/F domain
%           Df   background for each component to normalize the filtered raw data    

nA = sqrt(sum(A.^2))';
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = spdiags(nA,0,K,K)*C;

if nargin < 4
    [~,i] = min(sum(A.^6)); % identify background component
end

Yf = A'*(Y - A(:,[1:i-1,i+1:K])*C([1:i-1,i+1:K],:));
Df = median(Yf,2);
C_df = spdiags(Df,0,K,K)\C;
C_df(i,:) = 0;