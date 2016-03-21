function viewNeurons(obj, ind, C2, folder_nm)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components 
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

ind_del = false(size(ind));     % indicator of deleting neurons 
gSiz = obj.options.gSiz;        % maximum size of a neuron 

% time 
T = size(obj.C, 2); 
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

%% start viewing neurons 
for m=1:length(ind)
    %% full-frame view 
    subplot(221);
    imagesc(reshape(obj.A(:, ind(m)), obj.options.d1, obj.options.d2));
    axis equal; axis off;
    title(sprintf('Neuron %d', ind(m)));
    
    %% zoomed-in view 
    subplot(222);
    imagesc(reshape(obj.A(:, ind(m)), obj.options.d1, obj.options.d2));
    axis equal; axis off;
    x0 = center(ind(m), 2);
    y0 = center(ind(m), 1);
    xlim(x0+[-gSiz, gSiz]*2);
    ylim(y0+[-gSiz, gSiz]*2);
    
    
    %% temporal components 
    subplot(2,2,3:4);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(obj.A(:, ind(m)))); hold on;
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r');
    else
        
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
    end
    xlabel(str_xlabel); 
    
    %% save images 
    if save_img
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
    else
        fprintf('Neuron %d, type ''d'' to delete:   ', ind(m));
        temp = input('', 's');
        if temp=='d'
            ind_del(m) = true;
        end
    end
end
if save_img
    cd(cur_cd);
else
    obj.delete(ind(ind_del));
end

