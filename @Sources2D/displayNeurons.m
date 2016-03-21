function displayNeurons(obj, ind, C2, folder_nm)
%% view all components and delete components manually. it overlaps the 
%   neuon contours with correlation image and shows the corresponding 
%   spatial/temporal components 
%% input:  
%   ind: vector, indice of components to be displayed, no bigger than the maximum
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
center = obj.estCenter();       % estimate neurons center 
gSiz = obj.options.gSiz;        % maximum size of a neuron 
Cn = obj.Cn;                    % correlation image 
if isempty(Cn)
    fprintf('Please assign obj.Cn with correlation image!\n'); 
    return; 
end 
if isempty(obj.Coor)    % contours of the neuron has not been calculated 
    figure; obj.viewContours(obj.Cn, 0.8, 0); close; 
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
ctr = obj.estCenter(); 
for m=1:length(ind)
    %% contours + correlation image 
    subplot(221); cla;
    obj.image(Cn, [0,1]); hold on;
    axis equal off tight;
    % plot contour
    tmp_con = Coor{ind(m)};
    cont_del = (sum(tmp_con<=1, 1)>0); 
    tmp_con(:, cont_del) = [];
    if isempty(tmp_con)
        plot(ctr(m, 2), ctr(m, 2));
    else
        plot(tmp_con(1, 1:end), tmp_con(2, 1:end), 'r');
    end
    title(sprintf('Neuron %d', ind(m)));
    
    %% display spatial component 
    subplot(222);
    imagesc(reshape(obj.A(:, ind(m)), obj.options.d1, obj.options.d2));
    axis equal; axis off;
    x0 = center(ind(m), 2);
    y0 = center(ind(m), 1);
    xlim(x0+[-gSiz, gSiz]*2);
    ylim(y0+[-gSiz, gSiz]*2);
    
    %% display temporal traces 
    subplot(2,2,3:4);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(obj.A(:, ind(m)))); hold on;
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r');
    else
        
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
    end
    axis tight; 
    xlabel(str_xlabel); 
    
    %% save images 
    if save_img
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
    else
        fprintf('Neuron %d, type ''d'' to delete:   ', ind(m));
        temp = input('', 's');
        if temp=='d'
            ind_del(m) = true;
            disp(m); 
        end
    end
end
if save_img
    cd(cur_cd);
else
    obj.delete(ind(ind_del));
end

