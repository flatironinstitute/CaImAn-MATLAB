function str = cell2str(theCell)
% CELL2STR - Convert 1-D cells to a string
%
%   STR = CELL2STR(THECELL)
%
% Converts a 1-D cell to a string.
%
% Example: 
%   A = {'test','test2','test3'};
%   str = cell2str(A)
%
%       produced str = 
%    '{ 'test','test2','test3' }'
% 
 %1-dim cells only, only chars and matricies

if isempty(theCell), str = '{}'; return; end;

str = '{ ';
for i=1:length(theCell),
        if ischar(theCell{i})
                str = [str '''' theCell{i} ''', '];
        elseif isnumeric(theCell{i}),
                str = [str mat2str(theCell{i}) ', '];
        end;
end;
str = [str(1:end-2) ' }'];

