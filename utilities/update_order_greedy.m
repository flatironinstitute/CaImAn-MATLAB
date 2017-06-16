function [O,lo] = update_order_greedy(A)
    % greedy method for updating order
    K = size(A,2);
    if ~(size(A,1) == K)
        A = A'*A;
    end
    
    A = A>0;
    A(1:K+1:K^2) = 0;
    
    O{1} = 1;
    lo(1) = 1;
    for i = 2:K
       cont = true;
       for j = 1:length(O)           
           if ~any(A(i,O{j}))
               O{j} = [O{j},i];
               lo(j) = lo(j) + 1;
               cont = false;
               break;
           end           
       end
       if cont
           O{length(O)+1} = j;
           lo = [lo,1];
       end
    end    
end