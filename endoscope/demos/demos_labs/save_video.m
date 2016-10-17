Yac = neuron.reshape(neuron.A*neuron.C, 2);
Ybg = neuron.reshape(Ybg, 2);
% Y = neuron.reshape(Y, 2);
Y = neuron.reshape(Y, 2);


ctr = round( neuron.estCenter());
% figure;
% neuron.viewContours(Cn, .5, 0);
% cell_IDs = [];

figure('position', [0,0, 1024, 640]);
% avi_file = VideoWriter('~/Dropbox/Public/Garret/day1/residual.avi');
save_avi = true;
if save_avi
    avi_file = VideoWriter([dir_nm, file_nm, '_results1.avi']);
    avi_file.open();
end
Ymax = quantile(Y(1:1000:(d1*d2*T)), 0.999);
ACmax = quantile(Yac(1:1000:(d1*d2*T)), 0.9999);
for m=1:5:T
    subplot(4,6, [1,2, 5,6]);
    imagesc(Y(:, :,m), [0, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw data'); hold on; colorbar;
    
    subplot(4,6, [3,4, 7,8]);
    imagesc(Ybg(:, :, m), [0, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('background');
    colorbar;
    
        subplot(4,6, [9,10,13,14]);
        imagesc(Y(:, :, m)-Ybg(:, :, m), [0, ACmax]); hold on;
        set(gca, 'children', flipud(get(gca, 'children')));
        axis equal; axis off tight; title('raw-background'); hold on; colorbar;
    
        ax4 =  subplot(4,6, [11,12,15,16]);
        imagesc(neuron.reshape(neuron.A*neuron.C(:, m), 2), [0, ACmax]);
        %     imagesc(Ybg(:, :, m), [-50, 50]);
        text(d2/2-20, 10, sprintf('Time: %.2f Sec', m/neuron.Fs), 'color', 'w');
        axis equal off tight; colorbar; 
%     ax4 =  subplot(4,6, [11,12,15,16]);
%     imagesc(Y(:, :, m)-Ybg(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
%     set(gca, 'children', flipud(get(gca, 'children')));
%     axis equal; axis off tight; title('residual'); hold on; colorbar;
%     
%     subplot(4,6, [9,10,13,14]);
%     imagesc(neuron.reshape(neuron.A*neuron.C(:, m), 2), [0, ACmax]);
%     %     imagesc(Ybg(:, :, m), [-50, 50]);
%     text(d2/2-20, 10, sprintf('Time: %.2f Sec', m/neuron.Fs), 'color', 'w');
%     
    axis equal off tight; title('A\cdot C');
    colorbar;
    
    drawnow();
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [640, 1024]);
        avi_file.writeVideo(temp);
    end
end
if save_avi
    avi_file.close();
end

%% 
Yac = neuron.reshape(neuron.A*neuron.C, 2);
Ybg = neuron.reshape(Ybg, 2);
% Y = neuron.reshape(Y, 2);
Y = neuron.reshape(Y, 2);


ctr = round( neuron.estCenter());
% figure;
% neuron.viewContours(Cn, .5, 0);
% cell_IDs = [];

figure('position', [0,0, 1248, 600]);
% avi_file = VideoWriter('~/Dropbox/Public/Garret/day1/residual.avi');
save_avi = true;
if save_avi
    avi_file = VideoWriter([dir_nm, file_nm, '_results1.avi']);
    avi_file.open();
end
Ymax = quantile(Y(1:1000:(d1*d2*T)), 0.999);
ACmax = quantile(Yac(1:1000:(d1*d2*T)), 0.9999);

%     subplot(4,6, [5,6,11,12]); 
for m=1:5:T
    subplot(4,6, [1,2, 7, 8]);
    imagesc(Y(:, :,m), [0, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw data'); hold on; colorbar;
    
    subplot(4, 6, [1,2, 7, 8]+12);
         imagesc(Ybg(:, :, m), [0, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('background');
    colorbar;
    
    subplot(4,6, [3,4, 9,10]);
    imagesc(Y(:, :, m)-Ybg(:, :, m), [0, ACmax]); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw-background'); hold on; colorbar;
    
            
    subplot(4, 6, [3,4, 9,10]+12); 
    imagesc(Y(:, :, m)-Ybg(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('residual'); hold on; colorbar;

    subplot(4,6, [5,6,11,12]); 
        imagesc(neuron.reshape(neuron.A*neuron.C(:, m), 2), [0, ACmax]);
        %     imagesc(Ybg(:, :, m), [-50, 50]);
        text(d2/2-30, 10, sprintf('Time: %.2f Sec', m/neuron.Fs), 'color', 'w');
    axis equal off tight; title('A\cdot C'); colorbar; 
    
    drawnow();
    if save_avi
        temp = getframe(gcf);
        temp = imresize(temp.cdata, [600, 1248]);
        avi_file.writeVideo(temp);
    end
end
if save_avi
    avi_file.close();
end
%%
cell_id = [23, 95, 17, 110, 98 ];
nc = ceil((length(cell_id)+1)/2);
center = ctr(cell_id, :);
r0 = max(1, min(center(:, 1))-20); r1 = min(d1, max(center(:, 1))+20);
c0 = max(1, min(center(:, 2))-20); c1 = min(d2, max(center(:, 2))+20);
indr = r0:r1;
indc = c0:c1;
center = bsxfun(@minus, center, [r0, c0]-1);
Ysignal = neuron.reshape(Ysignal, 2); 
Ysignal_box = Ysignal(r0:r1, c0:c1, :);

figure('position', [100, 500, 1400, 480]);
avi_file = VideoWriter([dir_nm, file_nm, '_patch_neurons.avi']);
avi_file.open();
ACmax = 30;
ACmin = 5;
subplot(2, round(nc/2)+1, 1);
h_img = imagesc(Ysignal_box(:, :, 1), [ACmin, ACmax]); hold on;
axis equal off tight;
for m=1:length(cell_id)
    temp = text(center(m,2), center(m, 1), num2str(m), 'color', 'g', 'fontsize', 12, 'fontweight', 'bold');
    tmp_ctr = neuron.Coor{cell_id(m)};
    temp = plot(tmp_ctr(1, 2:end)-c0, tmp_ctr(2, 2:end)-r0, 'r');
end
for t=500:5:T
    subplot(2, nc, 1); hold off;
    imagesc(Ysignal_box(:, :, t), [ACmin, ACmax]);hold on;
    axis equal off tight;
    for m=1:length(cell_id)
        text(center(m,2), center(m, 1), num2str(m), 'color', 'g', 'fontsize', 12, 'fontweight', 'bold');
        tmp_ctr = neuron.Coor{cell_id(m)};
        plot(tmp_ctr(1, 5:end)-c0, tmp_ctr(2, 5:end)-r0, 'r');
    end
    title(sprintf('Time: %.2f Sec', t/neuron.Fs));
    
    for m=1:length(cell_id)
        subplot(2, nc, m+1);
        img = neuron.reshape(neuron.A(:, cell_id(m))*neuron.C(cell_id(m), t), 2);
        img = img(indr, indc);
        imagesc(img, [ACmin, ACmax]);
        axis equal off tight;
        title(sprintf('neuron %d', m));
    end
    drawnow();
    frame = getframe(gcf);
    frame.cdata = imresize(frame.cdata, [480, 1400]);
    avi_file.writeVideo(frame);
    % disp(t);
end
avi_file.close();





