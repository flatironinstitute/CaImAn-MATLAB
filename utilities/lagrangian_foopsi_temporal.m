function [c,ld] = lagrangian_foopsi_temporal(y,a,thr,G,ld_in)

% solves min sum(G*c) 
%       s.t.: G*c >= 0, ||a(j)*c - y(j)|| <= thr(j)
T = size(G,1);
la = length(a);
myfun = @(Ald) lagrangian_temporal_grad(Ald,y,thr);
options = optimset('GradObj','On','Display','Off','Algorithm','interior-point','TolX',1e-6);
if nargin < 5
    ld_in = 10*ones(length(a),1);
end
ld = fmincon(myfun,ld_in,[],[],[],[],zeros(length(a),1),[],[],options);

    function [f,grad] = lagrangian_temporal_grad(Al,y,thr)
        v = G'*ones(T,1);
        %options2 = optimset('Display','Off','Algorithm','interior-point-convex');
        %c = quadprog(2*sum(Al.*(a.^2))*speye(T),v-2*((la>1)*sum(spdiags(Al.*a,0,la,la)*y,1)' + (la==1)*y'*(Al.*a)),-G,zeros(T,1),[],[],[],[],[],options2);
        c = plain_foopsi(((la>1)*sum(spdiags(Al.*a,0,la,la)*y,1)' + (la==1)*y'*(Al.*a)-v/2)/(sum(Al.*(a.^2))),G);
        f = v'*c;    
        grad = (sum((a*c'-y).^2,2)-thr);
        f = f + Al(:)'*grad;
        f = -f;
        grad = -grad;
    end
end