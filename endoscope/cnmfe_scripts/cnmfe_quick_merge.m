neuron_bk = neuron.copy();
[merged_ROI, newIDs] = neuron.quickMerge(merge_thr);  % merge neurons based on the correlation computed with {'A', 'S', 'C'}
% A: spatial shapes; S: spike counts; C: calcium traces
if display_merge && ~isempty(merged_ROI)
    figure('position', [1,1, 1200, 600]);
    ind_before = false(size(neuron_bk.A, 2), 1);
    ind_after = false(size(neuron.A, 2), 1);
    m = 1;
    while m<=length(merged_ROI)
        subplot(221);
        tmp_img = neuron_bk.overlapA(merged_ROI{m});
        imagesc(tmp_img);
        axis equal off tight;
        subplot(222);
        imagesc(tmp_img);
        axis equal off tight;
        [tmp_r, tmp_c, ~] = find(sum(tmp_img, 3)>0);
        xlim([min(tmp_c)-10, max(tmp_c)+10]);
        ylim([min(tmp_r)-10, max(tmp_r)+10]);
        %         neuron.image(sum(Ain(:, merged_ROI{m}), 2));
        axis off;
        subplot(2,2,3:4);
        tmp_C = neuron_bk.C_raw(merged_ROI{m}, :)';
        tmp_C = bsxfun(@times, tmp_C, 1./max(tmp_C, [], 1));
        plot(tmp_C, 'linewidth', 2);
        
        temp = input('keep this merge? (y(default)/n(cancel)/b(back))/e(end)   ', 's');
        if strcmpi(temp, 'n')
            ind_after(newIDs(m)) = true;
            ind_before(merged_ROI{m}) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'e')
            break;
        else
            m = m+1;
        end
    end
    
    neuron.A = [neuron.A(:, ~ind_after), neuron_bk.A(:, ind_before)];
    neuron.C = [neuron.C(~ind_after, :); neuron_bk.C(ind_before, :)];
    neuron.C_raw = [neuron.C_raw(~ind_after, :); neuron_bk.C_raw(ind_before, :)];
    neuron.S = [neuron.S(~ind_after, :); neuron_bk.S(ind_before, :)];
    neuron.P.kernel_pars = [neuron.P.kernel_pars(~ind_after, :); neuron_bk.P.kernel_pars(ind_before, :)];
    clear neuron_bk;
end

% sort neurons
[Cpnr, srt] = sort(max(neuron.C, [], 2).*max(neuron.A, [], 1)', 'descend');
neuron.orderROIs(srt);
[Ain, Cin] = neuron.snapshot();   % keep the initialization results

%% view neurons
if view_neurons
    neuron.viewNeurons([], neuron.C_raw);
end