function [newSpikeTrain, newCalcium, newLL] = replaceSpike(oldSpikeTrain,oldCalcium,oldLL,filter,tau,obsCalcium,timetoRemove,indx,timetoAdd,Dt,A) %#codegen

% Replace a given spike with a new one in the existing spike train.

% Inputs:
% oldSpikeTrain:        current spike train
% oldCalcium:           current noiseless calcium trace
% oldLL:                current value of the log-likelihood function
% filters:              exponential rise and decay kernels for calcium transient
% tau:                  continuous time rise and decay time constants
% obsCalcium:           observed fluorescence trace
% timetoRemove:         time of the spike to be removed
% indx:                 place where the spike to be removed is in the existing spike train vector
% timetoAdd:            time of the spike to be added
% Dt:                   time-bin width
% A:                    spike amplitude

% Outputs:
% newSpikeTrain:        new vector of spike times
% newCalcium:           new noiseless calcium trace
% newLL:                new value of the log-likelihood function

% Author: Eftychios A. Pnevmatikakis and Josh Merel

    tau_h = tau(1);
    tau_d = tau(2);
    
    ef_h = filter{1,1};
    ef_d = filter{1,2};
    %ef_nh = filter{2,1};
    %ef_nd = filter{2,2};
    
    newSpikeTrain = oldSpikeTrain;
    newSpikeTrain(indx) = timetoAdd;
    
    %use infinite precision to scale the precomputed FIR approximation to the calcium transient    
    wk_hr = A*exp((timetoRemove - Dt*ceil(timetoRemove/Dt))/tau_h);
    wk_dr = A*exp((timetoRemove - Dt*ceil(timetoRemove/Dt))/tau_d);
    
    wk_ha = A*exp((timetoAdd - Dt*ceil(timetoAdd/Dt))/tau_h);
    wk_da = A*exp((timetoAdd - Dt*ceil(timetoAdd/Dt))/tau_d);

    %%%%%%%%%%%%%%%%%
    %handle ef_h first
    newCalcium = oldCalcium;
    min_t = floor(min(timetoRemove,timetoAdd));
    if any(ef_h)
        tmp = 1+ (min_t:min((length(ef_h)+min_t-1),length(newCalcium)-1));
        if floor(timetoRemove) == floor(timetoAdd)        
            wef_h = (wk_hr-wk_ha)*ef_h(1:length(tmp));
        elseif floor(timetoRemove) > floor(timetoAdd)
            wef_h = wk_hr*[zeros(1,min(length(tmp),floor(timetoRemove)-floor(timetoAdd))),ef_h(1:length(tmp)-(floor(timetoRemove)-floor(timetoAdd)))] - wk_ha*ef_h(1:length(tmp));
        else
            wef_h = wk_hr*ef_h(1:length(tmp)) - wk_ha*[zeros(1,min(length(tmp),floor(timetoAdd)-floor(timetoRemove))),ef_h(1:length(tmp)-(floor(timetoAdd)-floor(timetoRemove)))];
        end
        newCalcium(tmp) = newCalcium(tmp) - wef_h;
    end
    
    %%%%%%%%%%%%%%%%%
    %handle ef_d next
    tmp = 1+ (min_t:min((length(ef_d)+min_t-1),length(newCalcium)-1));
    if floor(timetoRemove) == floor(timetoAdd)        
        wef_d = (wk_dr-wk_da)*ef_d(1:length(tmp));
    elseif floor(timetoRemove) > floor(timetoAdd)
        wef_d = wk_dr*[zeros(1,min(length(tmp),floor(timetoRemove)-floor(timetoAdd))),ef_d(1:length(tmp)-(floor(timetoRemove)-floor(timetoAdd)))] - wk_da*ef_d(1:length(tmp));
    else
        wef_d = wk_dr*ef_d(1:length(tmp)) - wk_da*[zeros(1,min(length(tmp),floor(timetoAdd)-floor(timetoRemove))),ef_d(1:length(tmp)-(floor(timetoAdd)-floor(timetoRemove)))];
    end
    newCalcium(tmp) = newCalcium(tmp) - wef_d;

    obstemp = obsCalcium(tmp);
    newLL = oldLL - norm(newCalcium(tmp)-obstemp)^2 + norm(oldCalcium(tmp)-obstemp)^2; %+ 2*(newCalcium(tmp)-oldCalcium
    %%%%%%%%%%%%%%%%