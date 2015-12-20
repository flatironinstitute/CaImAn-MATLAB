function [c,b,c1,g,sn,sp] = MCEM_foopsi(y,b,c1,g,sn,options)
    defoptions.dt = 1;
    defoptions.MaxIter = 10;
    defoptions.MaxInerIter = 50;
    defoptions.TauStd = [0.2,2];
    defoptions.default_g = [0.6,0.9];
    
    if nargin < 6; options.defoptions; end 
    if ~isfield(options,'dt'); options.dt = defoptions.dt; end
    if ~isfield(options,'MaxIter'); options.MaxIter = defoptions.MaxIter; end
    if ~isfield(options,'MaxInerIter'); options.MaxInerIter = defoptions.MaxInerIter; end
    if ~isfield(options,'TauStd'); options.TauStd = defoptions.TauStd; end
    if ~isfield(options,'default_g'); options.default_g = defoptions.default_g; end
    
    dt = options.dt;
    
    if nargin < 5
        sn = [];
        if nargin < 4
            g = [];
            if nargin < 3
                c1 = [];
                if nargin < 2
                    b = [];
                end
            end
        end
    end

    
% initialization    
    
    [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,options);
    T = length(y);
    p = length(g);
    gr_ = sort(roots([1,-g(:)']));
    if p == 1; gr_ = [0,gr_]; end
    if any(gr_<0) || any(~isreal(gr_))
        gr_ = options.default_g;
    end
    tau = -dt./log(gr_);
    if p == 1        
        gr_(1) = 0;
        G1sp = zeros(T,1);
    else
        G1sp = make_G_matrix(T,gr_(1))\sp(:);
    end
    tau1_std = max(tau(1)/5,options.TauStd(1));
    tau2_std = min(tau(2)/10,options.TauStd(2));
    tauMoves = [0 0];
    tau_min = 0;
    tau_max = 2*tau(2);    
    gd_vec = max(gr_).^(0:T-1)';
    tau_sam = zeros(options.MaxIter,2);
    for iter = 1:options.MaxIter
        TAU = zeros(options.MaxInerIter,2);
        for iner_iter = 1:options.MaxInerIter
            if p >= 2       % update rise time constant
                logC = -norm(y(:) - c(:) - b - c1*gd_vec(:))^2; 
                tau_ = tau;
                tau_temp = tau_(1)+(tau1_std*randn); 
                while tau_temp >tau(2) || tau_temp<tau_min
                    tau_temp = tau_(1)+(tau1_std*randn);
                end 
                tau_(1) = tau_temp;
                gr_ = exp(dt*(-1./tau_));
                h_ = gr_(2)-gr_(1);
                G1_ = spdiags(ones(T,1)*[-min(gr_),1],-1:0,T,T);
                G1sp_ = G1_\sp(:);
                G2_ = spdiags(ones(T,1)*[-max(gr_),1],-1:0,T,T);
                c_ = (-G1sp_*gr_(1)+G2_\sp(:)*gr_(2))/h_;

                logC_ = -norm(y(:)-c_(:)-b-c1*gd_vec)^2;

                %accept or reject
                prior_ratio = 1;
        %         prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
                ratio = exp((logC_-logC)/(2*sn^2))*prior_ratio;
                if rand < ratio %accept
                    tau = tau_;
                    h = h_; G1sp = G1sp_; G2 = G2_; c = c_;
                    tauMoves = tauMoves + [1 1];
                else
                    tauMoves = tauMoves + [0 1];
                end                
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            % next update decay time constant
            %%%%%%%%%%%%%%%%%%%%%%%

            %initial logC
            logC = -norm(y(:)-c(:)-b-c1*gd_vec)^2;
            tau_ = tau;
            tau_temp = tau_(2)+(tau2_std*randn);
            while tau_temp>tau_max || tau_temp<tau_(1)
                tau_temp = tau_(2)+(tau2_std*randn);
            end  
            tau_(2) = tau_temp;

            gr_ = exp(dt*(-1./tau_));
            gd_vec_ = max(gr_).^(0:T-1)';
            h_ = gr_(2)-gr_(1);
            %G1_ = spdiags(ones(T,1)*[-min(gr_),1],[-1:0],T,T);
            G2_ = spdiags(ones(T,1)*[-max(gr_),1],[-1:0],T,T);  
            c_ = (-G1sp(:)*gr_(1)+G2_\sp(:)*gr_(2))/h_;

            logC_ = -norm(y(:)-c_(:)-b-c1*gd_vec_)^2;

            %accept or reject
            prior_ratio = 1;
    %       prior_ratio = gampdf(tau_(2),12,1)/gampdf(tau(2),12,1);
            ratio = exp((1./(2*sn^2)).*(logC_-logC))*prior_ratio;
            if rand<ratio %accept
                tau = tau_;
                %h = h_; G1 = G1_; G2 = G2_; 
                c = c_; gd_vec = gd_vec_;
                tauMoves = tauMoves + [1 1];
            else
                tauMoves = tauMoves + [0 1];
            end   
            TAU(iner_iter,:) = tau;
        end
        tau_ = mean(TAU);
        tau_sam(iter,:) = tau_;
        gr_ = exp(dt*(-1./tau_));
        %res = y(:)-c_(:)-b-c1*gd_vec;
        %sn   = 1./sqrt(gamrnd(1+T/2,1/(0.1 + sum((res.^2)/2))));
        g = [sum(gr_);-prod(gr_)];
        [c,b,c1,~,sn,sp] = constrained_foopsi(y,[],[],g(1:p),sn,options);
        %disp(tauMoves);
        
    end
    g.g = g;
    g.TAU = tau_sam(:,3-p:2);
end