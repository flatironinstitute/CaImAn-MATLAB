function [C,f,P,S] = update_temporal_components_parallel(Y,A,b,Cin,fin,P,options)

% update temporal components and background given spatial components in
% parallel, by forming of sequence of vertex covers. 
% Note that the parallel implementation is now implemented by default from
% update_temporal_components. update_temporal_components_parallel will be
% removed in the next release.

% Written by: 
% Eftychios A. Pnevmatikakis, Simons Foundation, 2015

if ~isempty(options); options.temporal_parallel = 1; end

[C,f,P,S] = update_temporal_components(Y,A,b,Cin,fin,P,options);