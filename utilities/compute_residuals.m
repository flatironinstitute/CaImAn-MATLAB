function R = compute_residuals(Y,A,b,C,f)

% compute residual traces for each component

% INPUTS
% Y     matrix of data in 2D format of pointer to file
% A     set of spatial components
% b     set of spatial background components
% C     set of temporal components
% f     set of temporal background component

AA = A'*A;
AY = mm_fun(A,Y);
nA2 = sum(A.^2,1);
R = bsxfun(@times, AY - AA*C - (A'*b)*f,1./nA2(:));