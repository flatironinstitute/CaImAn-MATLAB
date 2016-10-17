%%
Yac = neuron.reshape(neuron.A*neuron.C, 2);
Ybg = neuron.reshape(Ybg, 2);
% Y = neuron.reshape(Y, 2);
Ysignal = neuron.reshape(Ysignal, 2);
% ctr = round( neuron.estCenter());
% figure;
% neuron.viewContours(Cn, .5, 0);
% cell_IDs = [];

figure('position', [0,0, 1248, 600]);
% avi_file = VideoWriter('~/Dropbox/Public/Garret/day1/residual.avi');
if save_avi
    avi_file = VideoWriter([dir_nm, file_nm, '_results.avi']);
    avi_file.open();
end
temp  = quantile(Ybg(1:1000:numel(Ybg)), [0.0001, y_quantile]);
Ymin = temp(1);
Ymax = temp(2);
ACmax = quantile(Yac(1:1000:numel(Yac)), ac_quantile);

%     subplot(4,6, [5,6,11,12]);
for m=1:kt:T
    subplot(4,6, [1,2, 7, 8]);
    imagesc(Ybg(:, :,m)+Ysignal(:, :, m), [Ymin, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw data'); hold on; colorbar;
    
    subplot(4, 6, [1,2, 7, 8]+12);
    imagesc(Ybg(:, :, m), [Ymin, Ymax]);
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('background');
    colorbar;
    
    subplot(4,6, [3,4, 9,10]);
    imagesc(Ysignal(:, :, m), [0, ACmax]); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('raw-background'); hold on; colorbar;
    
    
    subplot(4, 6, [3,4, 9,10]+12);
    imagesc(Ysignal(:, :, m)-Yac(:, :, m), [-ACmax, ACmax]/2); hold on;
    set(gca, 'children', flipud(get(gca, 'children')));
    axis equal; axis off tight; title('residual'); hold on; colorbar;
    
    subplot(4,6, [5,6,11,12]);
    imagesc(Yac(:, :, m), [0, ACmax]);
    %     imagesc(Ybg(:, :, m), [-50, 50]);
    text(d2/2-30, 10, sprintf('Time: %.2f Sec', m/neuron.Fs), 'color', 'w');
    axis equal off tight; title('denoised'); colorbar;
    
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
