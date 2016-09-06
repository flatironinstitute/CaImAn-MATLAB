function ind_cell=run_movie(Y, A, C, Cn, min_max, Coor, ctr,  Ncell, tskip, save_avi, avi_name, S)
%% play movies together with selected neurons' calcium traces. It can be used for checking results
%% Input:
%
%% output

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
%% parameters
[d1,d2] = size(Cn);
if ismatrix(Y); Y=reshape(Y, d1, d2, []); T=size(Y, 3); else T=size(Y, ndims(Y)); end
if ~exist('save_avi', 'var');    save_avi=false;
elseif ~exist('avi_name', 'var')
    k=0; avi_name = sprintf('example_movie_%d', k);
    while exist(avi_name, file);  k=k+1; avi_name = sprintf('example_movie_%d', k); end
end
if ~exist('tskip', 'var')||isempty(tskip);    tskip=1; end
if ~exist('Ncell', 'var')||isempty(Ncell);    Ncell=5; end
if ~exist('Coor', 'var')||isempty(Coor); figure; Coor=plot_contours(A, Cn, 0.9, 0); close; end
if ~exist('ctr', 'var')|| isempty(ctr); ctr=com(A, d1, d2); end
if ~exist('min_max', 'var') || (isempty(min_max));
    temp = Y(:, :, randi(T, min(100, T), 1));
    min_max = quantile(temp(:), [0.2, 0.9999]);
    min_max(1) = max(min_max(1), 0);
end
ncol = 3;
%% select neurons by clicking the image.
figure('position', [100, 500, 1800, 420]);

% display correlation image for selecting neurons
subplot(Ncell, ncol, 1:ncol:(ncol*Ncell));
imagesc(Cn, [0, 1]);hold on;
axis equal; axis off;
ind_cell = zeros(Ncell,1);
m = 1;
while true
    subplot(Ncell, ncol, 1:ncol:(Ncell*ncol));
    
    [x, y] = ginput(1);
    h_dot = plot(x, y, '*m');
    dx = x-ctr(:, 2);
    dy = y-ctr(:, 1);
    dist_xy = sqrt(dx.^2+dy.^2);
    [~, ind_cell(m)] = min(dist_xy);
    temp = Coor{ind_cell(m)};
    h_con = plot(temp(1, 4:end), temp(2, 4:end), 'r');
    h_txt = text(ctr(ind_cell(m), 2), ctr(ind_cell(m), 1), num2str(m), 'fontsize', 10);
    disp(ctr(ind_cell(m), :)); 
    subplot(Ncell, ncol, ncol*(m-1)+[2, ncol]);cla;
    ylabel(num2str(m));
    plot(C(ind_cell(m), :)); hold on;
    ylabel(sprintf('ID %d', ind_cell(m))); 
    if (exist('S', 'var'))&&  (~isempty(S))
        tmp = find(S(ind_cell(m), :)>0);
        plot(tmp, S(ind_cell(m), tmp), '*g');
    end
    set(gca, 'xticklabel', []); axis tight;
    set(gca, 'yticklabel', []);
    xlim([1, T]);
    
    temp = input('keep the neuron? (y/n)   ', 's');
    if strcmp(temp, 'n')
        delete(h_dot); delete(h_con); delete(h_txt);
        continue;
    else
        ylabel(num2str(m));
        m =  m+1;
    end
    
    if m>Ncell; break; end
end
subplot(Ncell, ncol, ncol*Ncell+[-1,0]);
set(gca, 'xticklabel', get(gca, 'xtick'));
xlabel('Frame');

%% play movies
% create avi file if sve_avi==true
if save_avi
    avi_file = VideoWriter(avi_name);
    open(avi_file);
end
% draw selected neurons' contours and latex them
subplot(Ncell,ncol,1:ncol:(ncol*Ncell)); cla;  hold on; axis equal; axis off;
for m=1:Ncell
    temp = Coor{ind_cell(m)};
    plot(temp(1, 4:end), temp(2, 4:end), 'r');
    text(ctr(ind_cell(m), 2)+3, ctr(ind_cell(m), 1), num2str(m), 'fontsize', 10, 'color', 'g');
end
% play movie
for t=1:tskip:T
    subplot(Ncell,ncol,1:ncol:(ncol*Ncell));
    h_img = imagesc(Y(:, :, t), min_max);
%     colormap gray;
    hold on; axis equal; axis off;
    set(gca, 'children', flipud(get(gca, 'children')));
    
    % shift the red line indicating the current frame
    for m=1:Ncell
        subplot(Ncell, ncol, ncol*(m-1)+[2, ncol]);
        eval(sprintf('aa%d=plot([t,t], get(gca, ''ylim''), ''r''); ',m));
    end
    
    % save frame
    if save_avi
        frame = getframe(gcf);
        writeVideo(avi_file, frame);
    end
    drawnow();
    
    % remove red lines
    for m=1:Ncell
        eval(sprintf('delete(aa%d); ', m));
    end
    delete(h_img);
end
subplot(Ncell,ncol,1:ncol:(ncol*Ncell)); hold on; axis equal; axis off;
imagesc(Cn, [0, 1]);
set(gca, 'children', flipud(get(gca, 'children')));

if save_avi;  close(avi_file);end