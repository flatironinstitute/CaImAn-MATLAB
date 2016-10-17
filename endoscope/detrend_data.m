function [Ydt, X, R] = detrend_data(Y,nk)
%% detrend fluorescence signals with B-spline basis 
%% Inputs: 
%   Y: d X T matrix, video data 
%   nk: scalar, number of knots
% Outputs: 
%   Ydt: d X T matrix, detrended data 
%   X: T*M matrix, each column is one basis 
%   R: d*M matrix, coefficients of all basis for all pixels 

%% create basis 
[~, T] = size(Y); 

if ~exist('nk', 'var')
    nk = 5; 
end 
X = bsplineM((1:T)', linspace(1, T, nk), 4); 

%% compute coefficients of all spline basis 
R = (Y*X)/(X'*X); 

%% compute detrended data 
Ydt = Y-R*X'; 