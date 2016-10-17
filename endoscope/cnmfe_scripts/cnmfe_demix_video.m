%% play videos to verify the demixing results
cell_id = []; 
while true 
    temp = input('input cell ID (type 0 to end):     '); 
    temp = round(temp); 
    if and(temp>=1, temp<=size(neuron.A, 2))
        cell_id(end+1) = temp;  %#ok<SAGROW>
    else
        break; 
    end
end

nc = ceil((length(cell_id)+1)/2);
ctr = round( neuron.estCenter());
center = ctr(cell_id, :);
bd = 10; 
r0 = max(1, min(center(:, 1))-bd); r1 = min(d1, max(center(:, 1))+bd);
c0 = max(1, min(center(:, 2))-bd); c1 = min(d2, max(center(:, 2))+bd);
indr = r0:r1;
indc = c0:c1;
center = bsxfun(@minus, center, [r0, c0]-1);
Ysignal = neuron.reshape(Ysignal, 2);
Ysignal_box = Ysignal(r0:r1, c0:c1, :);

figure('position', [100, 500, 1400, 480]);
avi_file = VideoWriter([dir_nm, file_nm, '_patch_neurons_2.avi']);
avi_file.open();
ACmax = max(reshape(neuron.A(:, cell_id)*neuron.C(cell_id, :), 1, []))*0.5;
ACmin = 100;
subplot(2, round(nc/2)+1, 1);
h_img = imagesc(Ysignal_box(:, :, 1), [ACmin, ACmax]); hold on;
axis equal off tight;
for m=1:length(cell_id)
    text(center(m,2), center(m, 1), num2str(m), 'color', 'g', 'fontsize', 12, 'fontweight', 'bold');
    tmp_ctr = neuron.Coor{cell_id(m)};
    xx = tmp_ctr(1, :); yy = tmp_ctr(2, :); 
    ind = or(or(xx<c0, yy<r0), or(xx>c1, yy>r1)); 
    xx(ind) = []; 
    yy(ind) = []; 
    temp = plot(xx-c0, yy-r0, 'r');
end
for t=1:2:T
    subplot(2, nc, 1); hold off;
    imagesc(Ysignal_box(:, :, t), [ACmin, ACmax]);hold on;
    axis equal off tight;
    for m=1:length(cell_id)
    text(center(m,2), center(m, 1), num2str(m), 'color', 'g', 'fontsize', 12, 'fontweight', 'bold');
    tmp_ctr = neuron.Coor{cell_id(m)};
    xx = tmp_ctr(1, :); yy = tmp_ctr(2, :); 
    ind = or(or(xx<c0, yy<r0), or(xx>c1, yy>r1)); 
    xx(ind) = []; 
    yy(ind) = []; 
    temp = plot(xx-c0, yy-r0, 'r');
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
