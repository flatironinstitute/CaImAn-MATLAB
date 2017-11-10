function displayNeurons(obj, ind, C2, folder_nm)
%% view all components and delete components manually. it overlaps the
%   neuon contours with correlation image and shows the corresponding
%   spatial/temporal components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.
%   folder_nm: string, the folder to output images neuron by neuron.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save resulted images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end
%% delete neurons with too few pixels 
obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1)); 

ind_del = false(size(ind));     % indicator of deleting neurons
ctr = obj.estCenter();       % estimate neurons center
gSiz = obj.options.gSiz;        % maximum size of a neuron
Cn = obj.Cn;                    % correlation image
if isempty(Cn)
    fprintf('Please assign obj.Cn with correlation image!\n');
    return;
end
if isempty(obj.Coor) || (size(obj.A, 2)~=length(obj.Coor))   % contours of the neuron has not been calculated
    figure;
    obj.Coor = obj.get_contours();
end
Coor = obj.Coor;        % contours of all extracted neurons

% time
T = size(obj.C, 2);
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

% start displaying neurons
figure('position', [100, 100, 1024, 512]);
m = 1;
while and(m>=1, m<=length(ind))
    %% contours + correlation image
    subplot(221); cla; 
    obj.image(Cn, [0,1]); hold on; colormap winter;
    axis equal off tight;
    for k=1:m
        % plot contour
        tmp_con = Coor{ind(k)};
        cont_del = (sum(tmp_con<=1, 1)>0);
        tmp_con(:, cont_del) = [];
        if isempty(tmp_con)
            plot(ctr(m, 2), ctr(m, 2));
        else           
            if and(k<m, ~ind_del(k))
                plot(tmp_con(1, 1:end), tmp_con(2, 1:end), 'color','k', 'linewidth', 2);          
            elseif k==m
                plot(tmp_con(1, 1:end), tmp_con(2, 1:end), 'r', 'linewidth', 3);
            end
        end
        
        title(sprintf('Neuron %d', ind(m)));
    end
    
    %% zoomed-in view
    subplot(222);
    imagesc(reshape(obj.A(:, ind(m)), obj.options.d1, obj.options.d2));
    colormap(gca, jet);
    axis equal; axis off;
    x0 = ctr(ind(m), 2);
    y0 = ctr(ind(m), 1);
    xlim(x0+[-gSiz, gSiz]*2);
    ylim(y0+[-gSiz, gSiz]*2);
    
    
    %% temporal components
    subplot(2,2,3:4);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(obj.A(:, ind(m))), 'linewidth', 2); hold on;
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r');
    else
        
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
    end
    xlabel(str_xlabel);
    
    %% save images
    if save_img
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
    else
        fprintf('Neuron %d, keep(k, default)/delete(d)/split(s)/delete all(da)/backward(b)/end(e):    ', ind(m));
        tmp_option = input('', 's');
        if tmp_option=='d'
            ind_del(m) = true;
            m = m+1;
        elseif strcmpi(tmp_option, 'b')
            m = m-1;
        elseif strcmpi(tmp_option, 'da')
            ind_del(m:end) = true;
            break;
        elseif strcmpi(tmp_option, 'k')
            ind_del(m) = false;
            m= m+1;
        elseif strcmpi(tmp_option, 's')
            try
                subplot(222);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                tmpA = obj.A(:, ind(m));
                obj.A(:, end+1) = tmpA.*tmp_ind(:);
                obj.C(end+1, :) = obj.C(ind(m), :);
                obj.A(:, ind(m)) = tmpA.*(1-tmp_ind(:));
                obj.S(end+1, :) = obj.S(ind(m), :);
                obj.C_raw(end+1, :) = obj.C_raw(ind(m), :);
                obj.P.kernel_pars(end+1, :) = obj.P.kernel_pars(ind(m), :);
            catch
                sprintf('the neuron was not split\n');
            end
        elseif strcmpi(tmp_option, 't')
            try
                subplot(222);
                temp = imfreehand();
                tmp_ind = temp.createMask();
                obj.A(:, ind(m)) = obj.A(:, ind(m)).*tmp_ind(:);
            catch
                sprintf('the neuron was not trimmed\n');
            end
        elseif strcmpi(tmp_option, 'e')
            break;
        else
            m = m+1;
        end
    end
end
if save_img
    cd(cur_cd);
else
    obj.delete(ind(ind_del));
    obj.Coor = obj.get_contours(0.9); 
end
figure; obj.viewContours(obj.Cn, 0.8, 0); close;

