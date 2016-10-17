function lam = choose_lambda(kernel,sn,prob)
%% Choose regularizer weight for sparse deconvolution
%  The function chooses a regularization weight for sparse deconvolution
%  with a given kernel and noise level. The weight of the regularizer
%  is chosen such that noise alone cannot lead to a non-zero solution is
%  at least prob.

% Inputs:
% kernel:       deconvolution kernel 
%               if length(kernel) == 1 or length(kernel) == 2, then kernel
%               is treated as a set of discrete time constants g. Otherwise,
%               it is treated as the actual vector.
% sn:           noise level
% prob:         probability of zero solution (deafult: 0.99)

% Output:
% lam:          regularization weight

% Author: Eftychios A. Pnevmatikakis, 2016, Simons Foundation
% based on ideas and discussions with J. Tubiana and G. Debregeas, 
% Laboratorie Jean Parrin, UPMC, France    

% Reference for this approach:
% Selesnick, I. (2012). Sparse deconvolution (an MM algorithm)

%%

if nargin < 3 || isempty(prob)
    prob = 0.99;
end

if nargin < 2 || isempty(sn)
    sn = 1;
    warning('Noise value not provided. Using sn = 1...')
end

if length(kernel) <= 2
    kernel = filter(1,[1,-kernel(:)'],[1,zeros(1,999)]);
end

lam = sn*norm(kernel)*norminv(prob);