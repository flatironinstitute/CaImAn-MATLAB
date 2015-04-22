function [l c] = graph_connected_comp(sA)
% 
%  Computing connected components of an undirected graph - assuming sA is symmetric
% 
%  Usage:
%   [l c] = graph_conn_comp(sA);
% 
%  Inputs:
%   sA - sparse adjacency matrix (for directed graph - does not have to be symmetric)
% 
%  Outputs:
%   l - components labels
%   c - number of connected components
% 
% 
%  Compile using:
%  >> mex -O -largeArrayDims graph_conn_comp_mex.cpp
% 
% 

sA = spfun(@(x) ones(size(x)),sA);
if ~isequal(sA, sA')
    [ii jj] = find(sA);
    sA = sparse([ii jj],[jj ii], ones(1, 2*numel(ii)), size(sA,1), size(sA,2));
end
[l c] = graph_conn_comp_mex(sA);
l = double(l); % make it compatible of the rest of Matlab