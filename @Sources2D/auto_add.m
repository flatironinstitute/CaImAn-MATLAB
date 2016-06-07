function auto_add(obj, Yres, max_res)
%% manually add missing neurons

if nargin<3
    max_res = [];
end
if ~ismatrix(Yres); Yres = obj.reshape(Yres, 1); end;
Ymax = max(Yres, [], 2);
d1 = obj.options.d1;
d2 = obj.options.d2;
gSiz = obj.options.gSiz;
gSig = obj.options.gSig;

figure('position', [100, 500, 600, 400]);
subplot(221);
obj.image(Ymax);
axis equal off tight;

A = zeros(size(obj.A));
C = zeros(size(obj.C));
k = 0;  % number of neurons added
psf = fspecial('gaussian', round(gSiz), gSig);
psf = psf-mean(psf(:));
ind_search = true(d1*d2, 1);
while true
    % pick neuron locations
    subplot(221); cla;
    obj.image(Ymax.*ind_search); hold on;
    colorbar;
    axis equal off tight;
    [v_res, ind_ctr] = max(Ymax(:));
    if ~isempty(max_res) && v_res<max_res
        break;
    end
    [r, c] = ind2sub([d1,d2], ind_ctr);
    ind_search(ind_ctr) = false;
    plot(c, r, '*r');
    
    % select neighbours
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    Y_box = Yres(ind_nhood, :);
    HY_box = imfilter(reshape(Y_box, nr, nc,[]), psf, 'replicate');
    HY_box = reshape(HY_box, nr*nc, []);
    
    % run rank-1 NMF
    ind_ctr = sub2ind([nr, nc], r-rsub(1)+1, c-csub(1)+1);   % subscripts of the center
    sz = [nr, nc];
    [ai, ci, ind_success] =  extract_ac(HY_box, Y_box, ind_ctr, sz);
    if ~ind_success
        break;
    end
    subplot(222);
    imagesc(reshape(ai, nr, nc)); axis equal off tight;
    subplot(2,2,3:4);
    plot(ci);
    
    if ~isempty(max_res)
        k = k+1;
        A(ind_nhood, k) = ai;
        C(k, :) = ci;
        drawnow; 
    else
        temp = questdlg('save this neuron? ', '', 'Yes');
        if strcmpi(temp, 'Yes')
            k = k+1;
            A(ind_nhood, k) = ai;
            C(k, :) = ci;
        elseif strcmpi(temp, 'Cancel')
            break;
        end
    end
    Yres(ind_nhood, :) = Y_box - ai*ci;
    Ymax(ind_nhood) = max(Yres(ind_nhood, :), [], 2);
end

if k~=0
    obj.A = [obj.A, A(:, 1:k)];
    obj.C = [obj.C; C(1:k, :)];
end
























