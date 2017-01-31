function Y_ds = downsample_data(Y,direction,tsub,ssub,nrm)

% downsampling for 2d or 3d imaging data (already loaded in memory)

% INPUTS
% Y:            input dataset in 2D+T or 3D+T format
% direction:    direction of downsampling: 'space','time','spacetime' 
% tsub:         degree of downsampling in time
% ssub:         degree of downsampling in space
% nrm:          norm to be used when averaging (default: 1, plain averaging)

if nargin < 5 || isempty(nrm); nrm = 1; end
if nargin < 4 || isempty(ssub); ssub = 1; end
if nargin < 3 || isempty(tsub); tsub = 1; end
if nargin < 2 || isempty(direction); direction = 'spacetime'; end

ndimsY = ndims(Y) - 1;
if ndimsY == 2
    [d1,d2,T] = size(Y); d3 = 1;
elseif ndimsY == 3
    [d1,d2,d3,T] = size(Y);
end

d = [d1,d2,d3];
ds = d;

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
            if nrm == 1
                if ndimsY == 2; Y_ds = squeeze(nanmean(reshape(Y_ds(:, :, 1:(Ts*tsub)),ds(1), ds(2), tsub, Ts), 3)); end
                if ndimsY == 3; Y_ds = squeeze(nanmean(reshape(Y_ds(:, :, :, 1:(Ts*tsub)),ds(1), ds(2), ds(3), tsub, Ts), 4)); end
            elseif nrm == Inf
                if ndimsY == 2; Y_ds = squeeze(max(reshape(Y_ds(:, :, 1:(Ts*tsub)),ds(1), ds(2), tsub, Ts),[], 3)); end
                if ndimsY == 3; Y_ds = squeeze(max(reshape(Y_ds(:, :, :, 1:(Ts*tsub)),ds(1), ds(2), ds(3), tsub, Ts),[], 4)); end
            else
                if ndimsY == 2; Y_ds = squeeze(nanmean(reshape(Y_ds(:, :, 1:(Ts*tsub)).^nrm,ds(1), ds(2), tsub, Ts), 3)).^(1/nrm); end
                if ndimsY == 3; Y_ds = squeeze(nanmean(reshape(Y_ds(:, :, :, 1:(Ts*tsub)).^nrm,ds(1), ds(2), ds(3), tsub, Ts), 4)).^(1/nrm); end
            end
        end
    end
else
    fprintf('No downsampling is performed \n')
end