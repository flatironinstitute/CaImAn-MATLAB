function A_com = extract_patch(A,dims,patch_size,padding)

% Extractx a patch of size patch_size centered around the centroid of each
% component.
% INPUTS:
% A:            2d matrix of spatial components
% dims:         dimensions of FOV
% patch_size:   dimensions of patch
% padding:      if true components remain centered and are zero padded,
%               otherwise they are shifted (default: true)

% OUTPUT:
% A_com:            Nd matrix of patches

if ~exist('padding','var'); padding = true; end

nd = length(dims);
if nd == 2; dims(3) = 1; patch_size(3) = 1; end
K = size(A,2);
A = A/spdiags(sqrt(sum(A.^2,1))'+eps,0,K,K);      % normalize to sum 1 for each compoennt
cm = com(A,dims(1),dims(2),dims(3));
xx = -ceil(patch_size(1)/2-1):floor(patch_size(1)/2);
yy = -ceil(patch_size(2)/2-1):floor(patch_size(2)/2);
zz = -ceil(patch_size(3)/2-1):floor(patch_size(3)/2);
A_com = zeros([patch_size,K]);

for i = 1:K
    int_x = round(cm(i,1)) + xx;
    pad_pre_x = 0;
    pad_pre_y = 0;
    pad_post_x = 0;
    pad_post_y = 0;
    if int_x(1)<1
        if padding
            pad_pre_x = 1 - int_x(1);
            int_x = 1:int_x(end);
        else
            int_x= int_x + 1 - int_x(1);
        end
    end
    if int_x(end)>dims(1)
        if padding
            pad_post_x = int_x(end) - dims(1);
            int_x = int_x(1):dims(1);
        else
            int_x = int_x - (int_x(end)-dims(1));
        end
    end
    int_y = round(cm(i,2)) + yy;
    if int_y(1)<1
        if padding
            pad_pre_y = 1 - int_y(1);
            int_y = 1:int_y(end);
        else
            int_y = int_y + 1 - int_y(1);
        end
    end
    if int_y(end)>dims(2)
        if padding
            pad_post_y = int_y(end) - dims(2);
            int_y = int_y(1):dims(2);
        else
            int_y = int_y - (int_y(end)-dims(2));
        end
    end
    if nd == 3
        int_z = round(cm(i,3)) + zz;
        if int_z(1)<1
            int_z = int_z + 1 - int_z(1);
        end
        if int_z(end)>dims(3)
            int_z = int_z - (int_z(end)-dims(3));
        end
    else
        int_z = 1;
    end
    A_temp = reshape(full(A(:,i)),dims);
    A_temp = A_temp(int_x,int_y,int_z);
    if padding
        A_temp = padarray(A_temp,[pad_pre_x,pad_pre_y],0,'pre');
        A_temp = padarray(A_temp,[pad_post_x,pad_post_y],0,'post');
    end
    A_com(:,:,i) = A_temp;       
end