function AY = mm_fun(A,Y,chunk_size,run_paralel,FOV)
if ~exist('run_paralel', 'var') || isempty(run_paralel); run_paralel = true; end
if ~exist('chunk_size', 'var') || isempty(chunk_size); chunk_size = 2e4; end
if ~exist('FOV', 'var') || isempty(FOV); FOV = [512,512]; end

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

        if d1a == d1y    
            AY = zeros(d2a,d2y);
            if ~run_paralel
                for i = 1:chunk_size:d1a
                    AY = AY + A(i:min(i+chunk_size-1,d1a),:)'*double(Y.Yr(i:min(i+chunk_size-1,d1a),:));
                end
            else % adding parallel option
                chunks = 1:chunk_size:d1a; corenum = 4;
                nchunks = numel(chunks); % number of indexes
                chunk_idx = arrayfun(@(i) chunks(i:min((i+corenum-1), nchunks)), ...
                    1:corenum:nchunks, 'UniformOutput', false); % indices of the patches in each batch
                for i = 1:numel(chunk_idx)
                    batch2run = chunk_idx{i};
                    parfor ii = 1:numel(batch2run)
                        idx = batch2run(ii);
                        tempY = Y.Yr(idx:min(idx+chunk_size-1,d1a),:);
                        tempA = A(idx:min(idx+chunk_size-1,d1a),:);
                        AYt{ii} = tempA'*double(tempY);
                    end
                    AY = AY + sum(cat(3, AYt{:}),3);
                    clear AYt
                    if mod(i, 20) == 0
                        fprintf('%2.1f%% of chunks completed \n', i*100/numel(chunk_idx));
                    end
                end
            end
        elseif d1a == d2y
            AY = zeros(d1y,d2a);
            for i = 1:chunk_size:d2a
                AY(i:min(i+chunk_size-1,d1),:) = double(Y.Yr(i:min(i+chunk_size-1,d1a),:))*A;
            end
        elseif d2a == d1y
            AY = zeros(d1a,d2y);
            for i = 1:chunk_size:d2a
                AY = AY + A(:,i:min(i+chunk_size-1,d1a))'*double(Y.Yr(i:min(i+chunk_size-1,d1a),:));
            end
        elseif d2a == d2y
            AY = zeros(d1y,d1a);
            At = A';
            for i = 1:chunk_size:d1y
                AY(i:min(i+chunk_size-1,d1y),:) = double(Y.Yr(i:min(i+chunk_size-1,d1y),:))*At;
            end
        else
            error('matrix dimensions do not agree');
        end
    case 'hdf5'
        if d1a == d1y            
            AY = zeros(d2a,d2y);
            for t = 1:chunk_size:T
                Y_temp = bigread2(Y,t,min(T-t+1,chunk_size));
                Y_temp(isnan(Y_temp)) = 0;
                if ~ismatrix(Y_temp); Y_temp = reshape(Y_temp,d1y,[]); end
                AY(:,t:min(T,t+chunk_size-1)) = A'*double(Y_temp);
            end
        elseif d2a == d2y
            AY = zeros(d1y,d1a);            
            for t = 1:chunk_size:T
                Y_temp = bigread2(Y,t,min(T-t+1,chunk_size));
                Y_temp(isnan(Y_temp)) = 0;
                if ~ismatrix(Y_temp); Y_temp = reshape(Y_temp,d1y,[]); end
                AY = AY + double(Y_temp)*A(:,t:min(T,t+chunk_size-1))';
            end            
        end
end