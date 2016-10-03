function [Ain,Cin,bin,fin] = sparse_NMF_initialization(Y,K,options) %beta,eta,X0,err_thr,max_iter)

T = size(Y,ndims(Y));
if ~ismatrix(Y)
    Y = reshape(Y,numel(Y)/T,T);
end

defoptions = CNMFSetParms;
if nargin < 3 || isempty(options); options = defoptions; end
if ~isfield(options,'snmf_max_iter'); options.snmf_max_iter = defoptions.snmf_max_iter; end
    max_iter = options.snmf_max_iter;
if ~isfield(options,'err_thr'); options.err_thr = defoptions.err_thr; end
    err_thr = defoptions.err_thr;
if ~isfield(options,'eta'); options.eta = defoptions.eta; end
    eta = options.eta*max(Y(:))^2;
if ~isfield(options,'beta'); options.beta = defoptions.beta; end
    beta = options.beta;
if ~isfield(options,'nb'), options.nb = defoptions.nb; end 
    nb = options.nb;
if ~isfield(options,'rem_prct') || isempty(options.rem_prct); options.rem_prct = defoptions.rem_prct; end
    
C = rand(K,T);

repeat = 1;
iter = 1;
obj_ = 1e-10;

v = ver;
flag_optim = any(strcmp('Optimization Toolbox', {v.Name})); % check if optimization toolbox is present
if flag_optim
    if verLessThan('optim','6.3')
        min_options = optimset('Algorithm','interior-point','GradObj','On','Display','Off');
    else
        min_options = optimoptions('fmincon','Algorithm','interior-point','GradObj','On','Display','Off');
    end
end

% remove median
medY = prctile(Y,options.rem_prct,2);
Y = bsxfun(@minus, Y, medY);

while (iter <= max_iter) && repeat
    A = max((Y*C')*pinv(C*C'+beta*ones(K,K)),0);
    C = max((A'*A + eta*eye(K))\(A'*Y),0);
    
    ff = find(sum(C,2)==0);
    if ~isempty(ff)
        A(:,ff) = [];
        C(ff,:) = [];
        K = K - length(ff);
    end
    
    iter = iter + 1;
    if mod(iter,10) == 0;
        A = threshold_components(A,options);
        nC = double(sum(C.^2,2));
        AA = A'*A;
        mine = @(e) min_e(e,double(beta),double(eta),nC,A,AA);
        
        e = double(eta*nC./full(beta*A'*sum(A,2))).^(1/4);
        if flag_optim
            e = fmincon(mine,max(e,1e-4),[],[],[],[],1e-4*ones(K,1),[],[],min_options);
        end
        C = diag(e)\C;
        A = A*diag(e);
        fprintf('%i out of maximum %i iterations done \n',iter,max_iter);
    end
    obj = norm(Y - A*C,'fro')^2 + eta*norm(C,'fro')^2 + beta*norm(sum(A,2))^2;
    repeat = abs(obj - obj_) > err_thr*obj_;
    obj_ = obj;
end

fprintf('Algorithm converged after %i iterations. \n',iter-1);

[Ain,Cin] = order_components(A,C);
[bin,fin] = nnmf(max(Y - Ain*Cin + repmat(medY,1,T),0),nb);
Ain = sparse(double(Ain));
bin = double(bin);

    function [f,grad] = min_e(e,beta,eta,nC,A,AA)
        f = eta*norm(sqrt(nC)./e)^2 + beta*norm(A*e)^2;
        grad = -2*eta*nC./(e.^3) + 2*beta*AA*e;
    end


    function [A_or,C_or] = order_components(A,C)
        nA = sqrt(sum(A.^2));
        nr = length(nA);
        A = bsxfun(@times,A,1./nA(:)'); %A/spdiags(nA(:),0,nr,nr);
        mA = max(A);
        C = bsxfun(@times,C,nA(:)); %spdiags(nA(:),0,nr,nr)*C;
        nC2 = sqrt(sum(C.^2,2));
        mC = max(C,[],2);
        [~,srt] = sort(mC.*mA'./nC2,'descend');
        A_or = A(:,srt);
        C_or = C(srt,:);
    end

end