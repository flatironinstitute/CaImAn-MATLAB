function manual_add(obj, Yres)
%% manually add missing neurons

if ~ismatrix(Yres); Yres = obj.reshape(Yres, 1); end;
Ymax = max(Yres, [], 2);
d1 = obj.options.d1;
d2 = obj.options.d2;
gSiz = obj.options.gSiz;

figure;
subplot(221);
obj.image(Ymax);
axis equal off tight;

A = zeros(size(obj.A));
C = zeros(size(obj.C)); 
k = 0;  % number of neurons added 
while true
    % pick neuron locations
    subplot(221);
    obj.image(Ymax); 
    [c,r] = ginput(1); c=round(c); r=round(r);
    ind_ctr = sub2ind([d1,d2], r, c);
    y0 = squeeze(Yres(ind_ctr, :));
    
    % select neighbours
    rsub = max(1, -gSiz+r):min(d1, gSiz+r);
    csub = max(1, -gSiz+c):min(d2, gSiz+c);
    [cind, rind] = meshgrid(csub, rsub);
    [nr, nc] = size(cind);
    ind_nhood = sub2ind([d1, d2], rind(:), cind(:));
    Y_box = Yres(ind_nhood, :);
    
    % run rank-1 NMF
    ind_active = (corr(Y_box', y0')>0.8);
    [ai, ci] = nnmf(Y_box(ind_active,:), 1);
    subplot(222);
    temp = zeros(nr, nc);
    temp(ind_active) = ai;
    imagesc(temp); axis equal off tight;
    subplot(2,2,3:4);
    plot(ci);
    
    temp = questdlg('save this neuron? ', '', 'No'); 
    if strcmpi(temp, 'Yes')
        k = k+1; 
        A(ind_nhood(ind_active), k) = ai; 
        C(k, :) = ci; 
    elseif strcmpi(temp, 'Cancel')
        break; 
    end 
    Yres(ind_nhood(ind_active), :) = Y_box(ind_active, :) - ai*ci; 
    Ymax(ind_nhood) = max(Yres(ind_nhood, :), [], 2); 
end

if k~=0
    obj.A = [obj.A, A(:, 1:k)]; 
    obj.C = [obj.C; C(1:k, :)]; 
end 
























