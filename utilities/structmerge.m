function s_out = structmerge(s1, s2, varargin)
% STRUCTMERGE - Merge struct variables into a common struct
%
%  S_OUT = STRUCTMERGE(S1, S2, ...)
%
%  Merges the structures S1 and S2 into a common structure S_OUT
%  such that S_OUT has all of the fields of S1 and S2. When 
%  S1 and S2 share the same fieldname, the value of S2 is taken.
%
%  The behavior of the function can be altered by passing additional
%  arguments as name/value pairs. 
%
%  Parameter (default)     | Description
%  ------------------------------------------------------------
%  ErrorIfNewField (0)     | (0/1) Is it an error if S2 contains a
%                          |  field that is not present in S1?
% 
%  See also: STRUCT

ErrorIfNewField = 0;

assign(varargin{:});

f1 = fieldnames(s1);
f2 = fieldnames(s2);

if ErrorIfNewField,
	[c,f1i] = setdiff(f2,f1);
	if ~isempty(c),
		error(['Some fields of the second structure are not in the first: ' cell2str(c) '.']);
	end;
end;

s_out = s1;

for i=1:length(f2),
	s_out = setfield(s_out,f2{i},getfield(s2,f2{i}));
end;

