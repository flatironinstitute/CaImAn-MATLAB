function [W,H] = fast_nmf(A,B,k,max_iter)

% fast NMF computation of a matrix decomponsed in A*B form

if nargin < 4 || isempty(max_iter)
    max_iter = 20;
end

[W,H] = nnsvd(A,B,k);
H = H';

for iter = 1:max_iter
    if isempty(B)
        B = 1;
    end
    W = max((A*(B*H'))/(H*H'),0);
    H = max((W'*W)\((W'*A)*B),0);
    %disp(iter)
end

function [W,H] = nnsvd(A,B,k)

    mf = @(x,y) mat_vec(A,B,x,y);    
    sz = [size(A,1), size(B,2)];
    if isempty(B)
        sz(2) = size(A,2);
    end
    [U,S,V] = svds(mf,sz,k);

    W = zeros(sz(1),k);
    H = zeros(sz(2),k);

    for j = 1:k
        x = U(:,j);
        y = V(:,j);
        xp = (x>0).*x;
        xn = xp - x;
        yp = (y>0).*y;
        yn = yp - y;
        nxp = norm(xp);
        nxn = norm(xn);
        nyp = norm(yp);
        nyn = norm(yn);
        if nxp*nyp > nxn*nyn
            u = xp/nxp;
            v = yp/nyp;
            sigma = nxp*nyp;
        else
            u = xn/nxn;
            v = yn/nyn;
            sigma = nxn*nyn;
        end
        W(:,j) = sqrt(S(j,j)*sigma)*u;
        H(:,j) = sqrt(S(j,j)*sigma)*v;
    end
    if isempty(B)
        avg = mean(A(:));
    else
        avg = mean(A,1)*mean(B,2);
    end
    W(W==0) = avg;
    H(H==0) = avg;
    
    function y = mat_vec(A,B,x,method)

        if strcmp(method,'notransp')
            if isempty(B)
                y = A*x;
            else
                y = A*(B*x);
            end
        elseif strcmp(method,'transp')
            if isempty(B)
                y = A'*x;
            else
                y = B'*(A'*x);
                
            end
        end
    end
end

end