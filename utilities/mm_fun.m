function AY = mm_fun(A,Y,chunk_size)

% multiply A*Y or A'*Y or Y*A or Y*C' depending on the dimension for loaded
% or memory mapped Y.

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff');
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        sizY = size(Y);
        T = sizY(end);
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5');
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
    elseif strcmpi(ext,'raw')
        filetype = 'raw';
        fid = fopen(Y);
        FOV = [512,512];
        bitsize = 2;
        imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        T = file_length/imsize;
        sizY = [FOV,T];
        fclose(fid);
    end    
elseif isobject(Y);
    filetype = 'mem';
    sizY = size(Y,'Y');
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    Y = double(Y);
    if ~ismatrix(Y)
        Y = reshape(Y,[],size(Y,ndims(Y)));
    end
    sizY = size(Y);
    T = sizY(end);
end


[d1a,d2a] = size(A);

d1y = prod(sizY(1:end-1));
d2y = T;

switch filetype 
    case 'mat'
        if d1a == d1y   % A'*Y
            AY = A'*Y;
        elseif d1a == d2y
            AY = Y*A;
        elseif d2a == d1y
            AY = A*Y;
        elseif d2a == d2y  % Y*C'
            AY = Y*A';
        else
            error('matrix dimensions do not agree');
        end
    case 'mem'
        if nargin < 3 || isempty(chunk_size); chunk_size = 2e4; end
        if d1a == d1y
            if nargin < 3 || isempty(chunk_size); chunk_size = 2e4; end
            AY = zeros(d2a,d2y);
            for i = 1:chunk_size:d1a
                AY = AY + A(i:min(i+chunk_size-1,d1a),:)'*double(Y.Yr(i:min(i+chunk_size-1,d1a),:));
            end
        elseif d1a == d2y
            AY = zeros(d1y,d2a);
            for i = 1:chunk_size:d2a
                AY(i:min(i+chunk_size-1,d1),:) = double(Y.Yr(i:min(i+chunk_size-1,d1a),:))*A;
            end
        elseif d2a == d1y
            if nargin < 3 || isempty(chunk_size); chunk_size = 2e4; end
            AY = zeros(d1a,d2y);
            for i = 1:chunk_size:d2a
                AY = AY + A(:,i:min(i+chunk_size-1,d1a))'*double(Y.Yr(i:min(i+chunk_size-1,d1a),:));
            end
        elseif d2a == d2y
            AY = zeros(d1y,d1a);
            At = A';
            for i = 1:chunk_size:d1a
                AY(i:min(i+chunk_size-1,d1),:) = double(Y.Yr(i:min(i+chunk_size-1,d1a),:))*At;
            end
        else
            error('matrix dimensions do not agree');
        end
    case 'hdf5'
        if nargin < 3 || isempty(chunk_size); chunk_size = 2e3; end
        if d1a == d1y            
            AY = zeros(d2a,d2y);
            for t = 1:chunk_size:T
                Y_temp = bigread2(Y,t,min(T-t+1,chunk_size));
                if ~ismatrix(Y_temp); Y_temp = reshape(Y_temp,d1y,[]); end
                AY(:,t:min(T,t+chunk_size-1)) = A'*double(Y_temp);
            end
        elseif d2a == d2y
            AY = zeros(d1y,d1a);            
            for t = 1:chunk_size:T
                Y_temp = bigread2(Y,t,min(T-t+1,chunk_size));
                if ~ismatrix(Y_temp); Y_temp = reshape(Y_temp,d1y,[]); end
                AY = AY + double(Y_temp)*A(:,t:min(T,t+chunk_size-1))';
            end            
        end
end