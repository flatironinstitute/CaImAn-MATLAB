function [newSpikeTrain, newCalcium, newLL] = removeSpike(oldSpikeTrain,oldCalcium,oldLL,filters,tau,obsCalcium,timeToRemove,indx,Dt,A)

% Remove a given spike from the existing spike train.

% Inputs:
% oldSpikeTrain:        current spike train
% oldCalcium:           current noiseless calcium trace
% oldLL:                current value of the log-likelihood function
% filters:              exponential rise and decay kernels for calcium transient
% tau:                  continuous time rise and decay time constants
% obsCalcium:           observed fluorescence trace
% timetoRemove:         time of the spike to be removed
% indx:                 place where the spike to be removed is in the existing spike train vector
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
    
    newSpikeTrain = oldSpikeTrain;
    newSpikeTrain(indx) = [];
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient    
    wk_h = A*exp((timeToRemove - Dt*ceil(timeToRemove/Dt))/tau_h);
    wk_d = A*exp((timeToRemove - Dt*ceil(timeToRemove/Dt))/tau_d);
    
    %%%%%%%%%%%%%%%%%
    %handle ef_h first
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(ef_h)+floor(timeToRemove)-1),length(newCalcium)-1));
    wef_h = wk_h*ef_h(1:length(tmp));
    newCalcium(tmp) = newCalcium(tmp) - wef_h;

    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    relevantResidual(isnan(relevantResidual)) = 0;
    newLL = oldLL - ( wk_h^2*ef_nh(length(tmp)) + 2*relevantResidual*wef_h(:));
    oldCalcium = newCalcium;
    oldLL = newLL;
    %%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%
    %handle ef_d next
    newCalcium = oldCalcium;
    tmp = 1+ (floor(timeToRemove):min((length(ef_d)+floor(timeToRemove)-1),length(newCalcium)-1));
    wef_d = wk_d*ef_d(1:length(tmp));
    newCalcium(tmp) = newCalcium(tmp) - wef_d;

    relevantResidual = obsCalcium(tmp)-oldCalcium(tmp);
    relevantResidual(isnan(relevantResidual)) = 0;
    newLL = oldLL - ( wk_d^2*ef_nd(length(tmp)) + 2*relevantResidual*wef_d(:));
    %%%%%%%%%%%%%%%%