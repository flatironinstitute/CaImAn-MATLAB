function [Xs, bounce_count] = HMC_exact2(F, g, M, mu_r, cov, L, initial_X)

% Author: Ari Pakman

% returns samples from a d-dimensional Gaussian with m constraints given by  F*X+g >0 
% If cov == true
% then M is the covariance and the mean is mu = mu_r 
% if cov== false 
% then M is the precision matrix and the log-density is -1/2 X'*M*X + r'*X

% Input
% F:          m x d matrix
% g:          m x 1 vector 
% M           d x d matrix, must be symmmetric and definite positive
% mu_r        d x 1 vector. 
% cov:        see explanation above 
% L:          number of samples desired
% initial_X   d x 1 vector. Must satisfy the constraint.


% Output
% Xs:      d x L matrix, each column is a sample
% bounce_count:  number of times the particle bounced 

%% go to a whitened frame

m = size(g,1);
if size(F,1) ~= m
    display('error');
    return
end


if cov 
    mu= mu_r;
    g = g + F*mu;
    %if min(eig(M))<0
    %    M = M - 1.01*min(eig(M))*eye(size(M,1));
    %end
    R = chol(M);
    F = F*R';
    initial_X= initial_X -mu;
    initial_X = R'\initial_X;
    
else
    r=mu_r;
    R=chol(M);      % M = R'*R    
    mu = R\(R'\r);
    g = g+F*mu;
    F = F/R;               % this is the time-consuming step in this code section
    initial_X= initial_X -mu;
    initial_X = R*initial_X;    
end



    

d = size(initial_X,1);


Xs=NaN;

bounce_count =0;

nearzero= 10000*eps;




% Verify that initial_X is feasible
 c= F*initial_X +g;
 if any(c<0)
     display('error: inconsistent initial condition');
     return
 end


 
% Unsparcify the linear constraints
g=full(g);
F2 = sum(F.*F,2);  % squared norm of the rows of F, needed for reflecting the velocity
F=full(F);         % if we don't unsparcify  qj= F(j,:)*V/F2(j) becomes very slow.
Ft = F';
 

%% Sampling loop


last_X= initial_X;
Xs=zeros(d,L);
Xs(:,1)=initial_X;

i=2;
while (i <= L)
%i
stop=0;   
j=0;
V0= normrnd(0,1, d,1);   % initial velocity
X = last_X;

T=pi/2;                  % total time the particle will move
tt=0;                    % records how much time the particle already moved 

    while (1) 
        a = V0; 
        a= real(a);
        b = X;

        fa = F*a;
        fb = F*b;
        
        U = sqrt(fa.^2 + fb.^2);
        phi = atan2(-fa,fb);           % -pi < phi < +pi    



        pn = abs(g./U)<=1; % these are the walls that may be hit 


        % find the first time constraint becomes zero

        if any(pn) 
            inds = find(pn);

            phn= phi(pn);

            t1=-phn + acos(-g(pn)./U(pn));  % time at which coordinates hit the walls 
                                            % t1 in [-pi, 2*pi]
            t1(t1<0) = 2*pi + t1(t1<0);     % t1 in [0, 2*pi]                          
            t2 = -t1 -2*phn;                % second solution to hit the walls 
                                            % t2 in [-4*pi, 2*pi] 
            t2(t2<0) = 2*pi + t2(t2<0);     % t2 in [-2*pi, 2*pi] 
            t2(t2<0) = 2*pi + t2(t2<0);     % t2 in [0, 2*pi] 


         % if there was a previous reflection (j>0)
         % and there is a potential reflection at the sample plane                                    
         % make sure that a new reflection at j is not found because of numerical error
            if j>0    
                if pn(j)==1    
                    indj=sum(pn(1:j));            
                    tt1 = t1(indj);
                    if abs(tt1) < nearzero || abs(tt1-2*pi)< nearzero
                        t1(indj)=Inf;
                    else
                        tt2 = t2(indj);
                        if abs(tt2) < nearzero || abs(tt1-2*pi)< nearzero 
                            t2(indj) = Inf;
                        end
                    end                    
                end
            end


            [mt1 ind1] = min(t1);
            [mt2 ind2] = min(t2);

            [mt ind12]  = min([mt1 mt2]);

            if ind12==1
                m_ind = ind1;
            else
                m_ind= ind2;
            end

            % find the reflection plane 
            j = inds(m_ind);      % j is an index in the full vector of dim-m, not in the restriced vector determined by pn.

        else  %if pn(i) =0 for all i
                mt =T;
        end   %if pn(i)

        tt=tt+mt;
%        fprintf(num2str(tt/T));

        if tt>=T
            mt= mt-(tt-T);
            stop=1;

        end

        % move the particle a time mt

        X = a*sin(mt) + b*cos(mt);
        V = a*cos(mt) - b*sin(mt);

        if stop                    
            break;
        end

        % compute reflected velocity

        qj=  F(j,:)*V/F2(j);   
        V0 = V -2*qj*Ft(:,j);
        bounce_count = bounce_count+ 1;

        % dif =V0'*M*V0 - V'*M*V;

    end % while(1)

    % at this point we have a sampled value X, but due to possible
    % numerical instabilities we check that the candidate X satisfies the
    % constraints before accepting it.
    
    if all(F*X +g > 0)
        Xs(:,i)=X;
        last_X = X;
        i= i+1;
    
    else
    disp('hmc reject')
        
    end 

end %while (i <= L)

% transform back to the unwhitened frame
if cov
    Xs = R'*Xs  + repmat(mu, 1,L);   
else
    Xs = R\Xs  + repmat(mu, 1,L);   
end




end

