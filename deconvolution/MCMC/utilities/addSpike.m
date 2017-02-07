function [newSpikeTrain, newCalcium, newLL] = addSpike(oldSpikeTrain,oldCalcium,oldLL,filters,tau,obsCalcium,timeToAdd,indx,Dt,A)

% Add a given spike to the existing spike train.

% Inputs:
% oldSpikeTrain:        current spike train
% oldCalcium:           current noiseless calcium trace
% oldLL:                current value of the log-likelihood function
% filters:              exponential rise and decay kernels for calcium transient
% tau:                  continuous time rise and decay time constants
% obsCalcium:           observed fluorescence trace
% timetoAdd:            time of the spike to be added
% indx:                 place where the new spike is added in the existing spike train vector
% Dt:                   time-bin width
% A:                    spike amplitude

% Outputs:
% newSpikeTrain:        new vector of spike times
% newCalcium:           new noiseless calcium trace
% newLL:                new value of the log-likelihood function

% Author: Eftychios A. Pnevmatikakis and Josh Merel

    tau_h = tau(1);
    tau_d = tau(2);
    
    ef_h = filters{1,1};
    ef_d = filters{1,2};
    ef_nh = filters{2,1};
    ef_nd = filters{2,2};
    
    if isempty(oldSpikeTrain); indx = 1; end
    newSpikeTrain = [oldSpikeTrain(1:indx-1) timeToAdd oldSpikeTrain(indx:end)]; %possibly inefficient, change if problematic (only likely to be a problem for large numbers of spikes)
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient
    wk_h = A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau_h);
    wk_d = A*exp((timeToAdd - Dt*ceil(timeToAdd/Dt))/tau_d);
    
    %%%%%%%%%%%%%%%%%
    %handle ef_h first
    newCalcium = oldCalcium;
    tmp = 1 + (floor(timeToAdd):min((length(ef_h)+floor(timeToAdd)-1),length(newCalcium)-1));
    wef_h = wk_h*ef_h(1:length(tmp));
    newCalcium(tmp) = newCalcium(tmp) + wef_h;
    
    %if you really want to, ef*ef' could be precomputed and passed in
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    relevantResidual(isnan(relevantResidual)) = 0;
    %newLL = oldLL - ( wk_h^2*norm(ef_h(1:length(tmp)))^2 - 2*relevantResidual*(wk_h*ef_h(1:length(tmp))'));
    newLL = oldLL - ( wk_h^2*ef_nh(length(tmp)) - 2*relevantResidual*wef_h(:));
    oldCalcium = newCalcium;
    oldLL = newLL;
    %%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%
    %handle ef_d next
    tmp = 1 + (floor(timeToAdd):min((length(ef_d)+floor(timeToAdd)-1),length(newCalcium)-1));
    wef_d = wk_d*ef_d(1:length(tmp));
    newCalcium(tmp) = newCalcium(tmp) + wef_d;
    
    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    relevantResidual(isnan(relevantResidual)) = 0;
    %newLL = oldLL - ( wk_d^2*norm(ef_d(1:length(tmp)))^2 - 2*relevantResidual*(wk_d*ef_d(1:length(tmp))'));
    newLL = oldLL - ( wk_d^2*ef_nd(length(tmp)) - 2*relevantResidual*wef_d(:));
    %%%%%%%%%%%%%%%%%        