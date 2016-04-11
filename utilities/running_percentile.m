function [ y ] = running_percentile(x, win, p, varargin)
%RUNNING_PERCENTILE Median or percentile over a moving window.
%   Y = RUNNING_PERCENTILE(X,WIN,P) returns percentiles of the values in 
%   the vector Y over a moving window of size win, where p and win are
%   scalars. p = 0 gives the rolling minimum, and p = 100 the maximum.
%
%   running_percentile(X,WIN,P,THRESH) includes a threshold where NaN will be
%   returned in areas where the number of NaN values exceeds THRESH. If no
%   value is specified the default is the window size / 2.
%
%   The definition used is the same as for MATLAB's prctile: if the
%   window is too small for a given percentile, the max or min will be
%   returned. Otherwise linear interpolation is used on the sorted data.
%   Sorting is retained while updating the window, making this faster than
%   repeatedly calling prctile. The edges are handled by duplicating nearby
%   data.
%
%   Examples:
%      y = running_percentile(x,500,10); % 10th percentile of x with
%      a 500 point moving window

%   Author: Levi Golston, 2014

% Copyright (c) 2014, Levi Golston
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% Check inputs
N = length(x);
if win > N || win < 1
    error('window size must be <= size of X and >= 1')
end
if length(win) > 1
    error('win must be a scalar')
end
if p < 0 || p > 100
    error('percentile must be between 0 and 100')
end
if ceil(win) ~= floor(win)
    error('window size must be a whole number')
end
if ~isvector(x)
    error('x must be a vector')
end

if nargin == 4
    NaN_threshold = varargin{1};
else
    NaN_threshold = floor(win/2);
end

% pad edges with data and sort first window
if iscolumn(x)
    x = [x(ceil(win/2)-1 : -1 : 1); x; x(N : -1 : N-floor(win/2)+1); NaN];
else
    x = [x(ceil(win/2)-1 : -1 : 1), x, x(N : -1 : N-floor(win/2)+1), NaN];
end
tmp = sort(x(1:win));
y = NaN(N,1);

offset = length(ceil(win/2)-1 : -1 : 1) + floor(win/2);
numnans = sum(isnan(tmp));

% loop
for i = 1:N
	% Percentile levels given by equation: 100/win*((1:win) - 0.5);
	% solve for desired index
	pt = p*(win-numnans)/100 + 0.5;
    if numnans > NaN_threshold;   % do nothing
    elseif pt < 1        % not enough points: return min
		y(i) = tmp(1);
	elseif pt > win-numnans     % not enough points: return max
		y(i) = tmp(win - numnans);
	elseif floor(pt) == ceil(pt);  % exact match found
		y(i) = tmp(pt);
	else             % linear interpolation
		pt = floor(pt);
		x0 = 100*(pt-0.5)/(win - numnans);
		x1 = 100*(pt+0.5)/(win - numnans);
		xfactor = (p - x0)/(x1 - x0);
		y(i) = tmp(pt) + (tmp(pt+1) - tmp(pt))*xfactor;
    end
    
	% find index of oldest value in window
	if isnan(x(i))
		ix = win;  						  % NaN are stored at end of list
		numnans = numnans - 1;
	else
		ix = find(tmp == x(i),1,'first');
	end
	
	% replace with next item in data
	newnum = x(offset + i + 1);
	tmp(ix) = newnum;
	if isnan(newnum)
		numnans = numnans + 1;
	end
	tmp = sort(tmp);

end