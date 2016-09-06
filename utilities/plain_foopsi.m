function [Zin,ip_it] = plain_foopsi(H,D,I_est,eps)

% solves argmin ||X-H||^2 subject to D*X>=0 with an interior point method
% using I_est as the initial value and eps as the initial barrier weight


ln = length(H);
step_back_frac = 0.5;
iter = 0;
if nargin == 2
    I_est = 1e-3*ones(ln,1);
    eps = 1;
end
Zin = I_est(:);

if nargin == 3
    eps = 1;
end
while eps>1e-5
    n = D*Zin;
    nnd = 10;
    E = norm(Zin-H)^2 - eps*sum(log(D*Zin));
    grad = 2*(Zin-H) - eps*D'*(n.^(-1));
    Hs = 2*speye(ln) + eps*D'*spdiags(n.^(-2),0,ln,ln)*D;          
    while nnd/2>1
        iter = iter + 1;
        Z_dir = -Hs\grad;
        hit = -n./(D*Z_dir);
        if all(hit<0)
            s = 1;
        else
            s = min(1,.9*min(hit(hit>=0)));
        end
        E_new = E; s = s/step_back_frac;
        x_dg = grad'*Z_dir;
        while E_new > E + 0.25*s*x_dg
            s=s*step_back_frac; 
            Z_new = Zin + s*Z_dir;
            n = D*Zin;
            E_new = norm(Z_new-H)^2 - eps*sum(log(D*Z_new));
        end
        %E = E_new;
        Zin = Zin + s*Z_dir;
        nnd = -x_dg;
        E = norm(Zin-H)^2 - eps*sum(log(D*Zin));
        n = D*Zin;
        grad = 2*(Zin-H) - eps*D'*(n.^(-1));
        Hs = 2*speye(ln) + eps*D'*spdiags(n.^(-2),0,ln,ln)*D;
        %disp(nnd)
    end
    eps = eps/10;
end
%fprintf('Interior point method converged after %i iterations \n',iter);
ip_it = iter;
