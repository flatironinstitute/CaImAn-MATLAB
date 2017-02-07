function [new_spikes, new_calcium, moves]  = get_next_spikes(curr_spikes,curr_calcium,calciumSignal,filters,tau,calciumNoiseVar, lam, proposalSTD, add_move, Dt, A, con_lam)

% Function for sampling the next of spike times given the current set of spike time and observed fluorescence trace

% Inputs:
% curr_spikes:          current set of spike times in continuous time
% curr_calcium:         current estimate of noiseless calcium trace
% calciumSingal:        observed fluorescence trace
% filters:              exponential rise and decay kernels for calcium transient
% tau:                  continuous time rise and decay time constants
% calciumNoiseVar:      observation noise variance
% lam:                  function handle for computing firing rate
% proposalSTD:          standard deviation for perturbing a spike time          
% add_move:             # of addition/removal proposals per sample
% Dt:                   time-bin width
% A:                    spike amplitude
% con_lam:              flag for constant firing rate (speeds up computation)

% Outputs:
% new_spikes:           new set of spike times
% new_calcium:          new estimate of noiseless calcium trace
% moves:                acceptance probabilities for perturbing,adding, droping spikes respectively

% Author: Eftychios A. Pnevmatikakis and Josh Merel

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% initialize some parameters
    T = length(calciumSignal); %for all of this, units are bins and spiketrains go from 0 to T where T is number of bins
    ff = ~isnan(calciumSignal); % entries with data
       
    %% start with initial spiketrain and initial predicted calcium 
    si = curr_spikes; %initial set of spike times has no spikes - this will not be sorted but that shouldn't be a problem
    new_calcium = curr_calcium; %initial calcium is set to baseline 
    
    N = length(si); %number of spikes in spiketrain
        
    %initial logC - compute likelihood initially completely - updates to likelihood will be local
    logC = -norm(new_calcium(ff)-calciumSignal(ff))^2; 
       
    %flag for uniform vs likelihood proposal (if using likelihood proposal, then time shifts are pure Gibbs)
    %this really should be split into four cases (not implemented yet), 
    % 1) RW for time shifts with uniform add/drop
    % 2) RW for time shifts with likelihood proposal add/drop
    % 3) Gibbs for time shifts with uniform add/drop
    % 4) Gibbs for time shifts with likelihood proposal add/drop
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    addMoves = [0 0]; %first elem is number successful, second is number total
    dropMoves = [0 0];
    timeMoves = [0 0];
    
    %% loop over spikes, perform spike time move (could parallelize here for non-interacting spikes, i.e. spikes that are far enough away)
        
    for ni = 1:N %possibly go through in a random order (if you like)

        tmpi = si(ni);
        tmpi_ = si(ni)+(proposalSTD*randn); %with bouncing off edges
        if tmpi_<0
            tmpi_ = -(tmpi_);
        elseif tmpi_>T
            tmpi_ = T-(tmpi_-T);
        end           

        %set si_ to set of spikes with the move and ci_ to adjusted calcium and update logC_ to adjusted
%           [si_, ci_, logC_] = removeSpike(si,ci,logC,ef,tau,calciumSignal,tmpi,ni,Dt,A);
%           [si_, ci_, logC_] = addSpike(si_,ci_,logC_,ef,tau,calciumSignal,tmpi_,ni,Dt,A);     
        [si_, ci_, logC_] = replaceSpike(si,new_calcium,logC,filters,tau,calciumSignal,tmpi,ni,tmpi_,Dt,A);

        %accept or reject
        ratio = exp((logC_-logC)/(2*calciumNoiseVar));
        if ~con_lam; ratio = ratio*lam(tmpi)/lam(tmpi_); end
        if ratio>1 %accept
            si = si_;
            new_calcium = ci_;
            logC = logC_;
            timeMoves = timeMoves + [1 1];
        elseif rand<ratio %accept
            si = si_;
            new_calcium = ci_;
            logC = logC_;
            timeMoves = timeMoves + [1 1];
        else
            %reject - do nothing
            timeMoves = timeMoves + [0 1];
        end

    end
    N = length(si);
                       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% loop over add/drop a few times
    %define insertion proposal distribution as the likelihood function
    %define removal proposal distribution as uniform over spikes
    %perhaps better is to choose smarter removals.
    for ii = 1:add_move 
        %% add
        %propose a uniform add
        tmpi = T*Dt*rand;         
        [si_, ci_, logC_] = addSpike(si,new_calcium,logC,filters,tau,calciumSignal,tmpi,length(si)+1,Dt,A);

        %forward probability
        fprob = 1/(T*Dt);

        %reverse (remove at that spot) probability
        rprob = 1/(N+1);

        %accept or reject
        ratio = exp((1/(2*calciumNoiseVar))*(logC_ - logC))*(rprob/fprob)*lam(tmpi); %posterior times reverse prob/forward prob
        if ratio>1 %accept
            si = si_;
            new_calcium = ci_;
            logC = logC_;
            addMoves = addMoves + [1 1];
        elseif rand<ratio %accept
            si = si_;
            new_calcium = ci_;
            logC = logC_;
            addMoves = addMoves + [1 1];
        else
            %reject - do nothing
            addMoves = addMoves + [0 1];
        end
        N = length(si);

        %% delete
        if N>0                
            %propose a uniform removal
            tmpi = randi(N);
            [si_, ci_, logC_] = removeSpike(si,new_calcium,logC,filters,tau,calciumSignal,si(tmpi),tmpi,Dt,A);

            %reverse probability
            rprob = 1/(T*Dt);

            %compute forward prob
            fprob = 1/N;

            %accept or reject
            ratio = exp((logC_ - logC)/(2*calciumNoiseVar))*(rprob/fprob)*(1/lam(si(tmpi))); %posterior times reverse prob/forward prob

            if ratio>1 %accept
                si = si_;
                new_calcium = ci_;
                logC = logC_;
                dropMoves = dropMoves + [1 1];
            elseif rand<ratio %accept
                si = si_;
                new_calcium = ci_;
                logC = logC_;
                dropMoves = dropMoves + [1 1];
            else
                %reject - do nothing
                dropMoves = dropMoves + [0 1];
            end
            N = length(si);

        end                
    end
    new_spikes = si;
    moves = [timeMoves; addMoves; dropMoves];    