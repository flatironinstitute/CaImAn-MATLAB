function savefast(filename, varargin)
% savefast: fast saves of large arrays to .mat files
%
% Matlab's 'save' command can be very slow when saving large arrays,
% because by default Matlab attempts to use compression. This function
% provides a much faster alternative, at the cost of larger files.
%
% The syntax is identical to that of the Matlab save command.
%
% Example:
% >> ops = struct('algorithm', 'greedy');
% >> A = int32(randi(20, 1000, 1200, 40));
% >> B = randn(500, 1800, 60);
% >> tic; save /tmp/test ops A B; toc
% Elapsed time is 22.980294 seconds.
% >> tic; savefast /tmp/test ops A B; toc
% Elapsed time is 0.571098 seconds.

% Copyright 2013 by Timothy E. Holy
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

  % Extract the variable values
  vars = cell(size(varargin));
  for i = 1:numel(vars)
    vars{i} = evalin('caller', varargin{i});
  end
  
  % Separate numeric arrays from the rest
  isnum = cellfun(@(x) isa(x, 'numeric'), vars);
  
  % Append .mat if necessary
  [filepath, filebase, ext] = fileparts(filename);
  if isempty(ext)
    filename = fullfile(filepath, [filebase '.mat']);
  end
  
  create_dummy = false;
  if all(isnum)
    % Save a dummy variable, just to create the file
    dummy = 0; %#ok<NASGU>
    save(filename, '-v7.3', 'dummy');
    create_dummy = true;
  else
    s = struct;
    for i = 1:numel(isnum)
      if ~isnum(i)
        s.(varargin{i}) = vars{i};
      end
    end
    save(filename, '-v7.3', '-struct', 's');
  end
  
  % Delete the dummy, if necessary, just in case the user supplied a
  % variable called dummy
  if create_dummy
    fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    H5L.delete(fid,'dummy','H5P_DEFAULT');
    H5F.close(fid);
  end
  
  % Save all numeric variables
  for i = 1:numel(isnum)
    if ~isnum(i)
      continue
    end
    varname = ['/' varargin{i}];
    h5create(filename, varname, size(vars{i}), 'DataType', class(vars{i}));
    h5write(filename, varname, vars{i});
  end
end

