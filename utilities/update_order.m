function [O,lo] = update_order(A)

    K = size(A,2);
    F = (A'*A>0);       % find overlapping components
    F(1:K+1:K^2) = 0;   % remove diagonal elements
    rem_ind = 1:K;      % remaining indeces

    dp = 0;
    while ~isempty(rem_ind);
        dp = dp+1;
        L = sort(app_vertex_cover(F(rem_ind,rem_ind)),'ascend');
        ord_ind = setdiff(rem_ind,rem_ind(L));
        O{dp} = ord_ind;
        lo(dp) = length(ord_ind);
        rem_ind = rem_ind(L);
    end
    O = fliplr(O);
    lo = fliplr(lo);
    
    function L = app_vertex_cover(A)
        L = [];
        while ~isempty(find(A, 1))
            i = randsample(find(A),1);
            [u,~] = ind2sub(size(A),i);
            L = [L,u];
            A(u,:) = 0; % edges from u
            A(:,u) = 0; % edges to u
        end
    end
end