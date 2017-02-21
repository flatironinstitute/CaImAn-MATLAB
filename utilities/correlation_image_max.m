function CI = correlation_image_max(Y,sz,sizY,batch_size,min_batch_size)

% computes the correlation image in batches and then taking the max

% Y:                dataset loaded in memory or memory mapped file
% sz:               size of neighborhood
% sizY:             spatial dimensions of dataset
% batch_size:       size of each batch
% min_batch_size:   minimum allowed batch_size

memmaped = isobject(Y);

if ~exist('min_batch_size', 'var') || isempty(min_batch_size); min_batch_size = 1000; end
if ~exist('batch_size', 'var') || isempty(batch_size); batch_size = 2000; end
min_batch_size = min(min_batch_size,batch_size);
if ~exist('sizY','var') || isempty(sizY); 
    if memmaped; sizY = Y.sizY; else sizY = size(Y); end 
     T = sizY(end); sizY = sizY(1:end-1);
else
    if memmaped; T = Y.sizY(end); else T = size(Y,ndims(Y)); end
end
nd = length(sizY);

if ~exist('sz','var'); 
    if nd == 2;  sz = 4; elseif nd == 3; sz = 6; end
end

ln = max(floor(T/batch_size) + (rem(T,batch_size)>min_batch_size),1);

Cn = cell(ln,1);

for i = 1:ln-1
    if nd == 2
        if memmaped
            Cn{i} = correlation_image_3D(single(Y.Y(:,:,(i-1)*batch_size+1:i*batch_size)),sz,sizY);
        else
            Cn{i} = correlation_image_3D(single(Y(:,:,(i-1)*batch_size+1:i*batch_size)),sz,sizY);
        end
    elseif nd ==3 
        if memmaped
            Cn{i} = correlation_image_3D(single(Y.Y(:,:,:,(i-1)*batch_size+1:i*batch_size)),sz,sizY);
        else
            Cn{i} = correlation_image_3D(single(Y(:,:,:,(i-1)*batch_size+1:i*batch_size)),sz,sizY);
        end
    end
end

if nd == 2
    if memmaped
        Cn{ln} = correlation_image_3D(single(Y.Y(:,:,(ln-1)*batch_size+1:T)),sz,sizY);
    else
        Cn{ln} = correlation_image_3D(single(Y(:,:,(ln-1)*batch_size+1:T)),sz,sizY);
    end
elseif nd ==3 
    if memmaped
        Cn{ln} = correlation_image_3D(single(Y.Y(:,:,:,(ln-1)*batch_size+1:T)),sz,sizY);
    else
        Cn{ln} = correlation_image_3D(single(Y(:,:,:,(ln-1)*batch_size+1:T)),sz,sizY);
    end
end

CI = max(cat(nd+1,Cn{:}),[],nd+1);