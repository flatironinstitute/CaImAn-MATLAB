function idx = patch_to_linear(patch, sizY)
    % helper function to build linear indices from patch start/stop indices
    slice_idx = patch_to_indices(patch);
    subs_idx = cell(1, numel(slice_idx));
    [subs_idx{:}] = ndgrid(slice_idx{:});
    subs_idx = cellfun(@(x) x(:), subs_idx, 'un', false);
    idx = sub2ind(sizY(1:end-1), subs_idx{:});
    
    function idx = patch_to_indices(patch)
        % helper function to build indices vector from patch start/stop indices
        idx = arrayfun(@(x,y) x:y, patch(1:2:end), patch(2:2:end), 'un', false);
    end
end