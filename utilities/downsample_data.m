function [Y_ds,options_ds] = downsample_data(Y,direction,options)

% downsampling for 2d imaging data

defoptions = CNMFSetParms;
if nargin < 3 || isempty(options);
    options = defoptions;
end

if nargin < 2 || isempty(direction)
    direction = 'spacetime';
end

if ismatrix(Y)
    if options.d3 == 1;
        Y = reshape(Y,options.d1,options.d2,[]);
    else
        Y = reshape(Y,options.d1,options.d2,options.d3,[]);
    end
end

if ~isfield(options, 'ssub'); options.ssub = defoptions.ssub; end; ssub = options.ssub;
if ~isfield(options, 'tsub'), options.tsub = defoptions.tsub; end; tsub = options.tsub;

d = [options.d1,options.d2,options.d3];
ds = d;
ndimsY = ndims(Y)-1;
T = size(Y,ndimsY+1);
Ts = floor(T/tsub);          %reduced number of frames
Y_ds = Y;
if ~((ssub == 1) && (tsub == 1))
    fprintf('starting resampling \n')
    if strcmpi(direction,'space') || strcmpi(direction,'spacetime');
        if ssub~=1;
            ds(1:2) = ceil(d(1:2)/ssub); % do not subsample along z axis
            if ndimsY == 2; Y_ds = imresize(Y_ds, [ds(1), ds(2)], 'box'); end
            if ndimsY == 3;
                Y_ds = zeros([ds(1:2),T,ds(end)]);
                for z = 1:ds(3)
                    Y_ds(:,:,:,z) = imresize(squeeze(Y_ds(:,:,z,:)), [ds(1), ds(2)], 'box');
                end
                Y_ds = permute(Y_ds,[1,2,4,3]);
            end
        else
            Y_ds = Y;
        end
    end
    if strcmpi(direction,'time') || strcmpi(direction,'spacetime')
        if tsub~=1
            if ndimsY == 2; Y_ds = squeeze(mean(reshape(Y_ds(:, :, 1:(Ts*tsub)),ds(1), ds(2), tsub, Ts), 3)); end
            if ndimsY == 3; Y_ds = squeeze(mean(reshape(Y_ds(:, :, :, 1:(Ts*tsub)),ds(1), ds(2), ds(3), tsub, Ts), 4)); end
        end
    end
else
    fprintf('No downsampling is performed \n')
end

options_ds = options;
options_ds.d1 = ds(1);
options_ds.d2 = ds(2);