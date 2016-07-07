function [ inferred, filtered, raw ] = signalExtraction(Y,A,C,b,f,d1,d2,extractControl)
% this code extract the signal after CNMF is ran 

% inputs: Y raw data (d X T matrix, d # number of pixels, T # of timesteps)
%         A matrix of spatial components (d x K matrix, K # of components)
%         C matrix of temporal components (K x T matrix)
%         S matrix of deconvolved activity ((K-1) x T matrix) (optional)
%         b matrix of background spatial components (d x options.nb, options.nb # of components)
%         f matrix of background temporal components (options.nb x T matrix, options.nb # of components) 
%         extractControl.baselineRatio: residue baseline of signals from CNMF (default:0.25); set to 0 to ignore this residue baseline 
%         extractControl.thr: energy fraction included in the raw signal contour (default:0.99)

% inferred: spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
% filtered: spatial weighting on the ROI pixels, unmixing, background substraction
% raw: uniform spatial weighting on the ROI pixels (with threshold to remove the very low energy pixels), background substraction
% df: df signal
% F: baseline
% dfof: dfof
% CR: information of the cell contour used in raw signal calculation

% Author: Eftychios A. Pnevmatikakis, and Weijian Yang

if isfield(extractControl,'baselineRatio') baselineRatio=extractControl.baselineRatio;
else baselineRatio=0.25; end

if isfield(extractControl,'thr') thr=extractControl.thr;
else thr=0.99; end

% normalize spatial components to unit energy
A2 = [A,b];
C2 = [C;f];
K = size(C,1);
K2 = size(C2,1);
nA = sqrt(sum(A2.^2))';
A2 = A2/spdiags(nA,0,K2,K2);    % normalize spatial components to unit energy
C2 = spdiags(nA,0,K2,K2)*C2;
A = A2(:,1:K);
C = C2(1:K,:);
b = A2(:,K+1:end);
f = C2(K+1:end,:);

% infer signal
Yf = A'*(Y - A*C);
residueF = getBaseline(C,baselineRatio);
inferred.F = median(Yf,2) + residueF;
inferred.df = C - repmat(residueF,1,size(C,2));
inferred.dfof = spdiags(inferred.F,0,K,K)\inferred.df;

% filtered signal
filtered.F = inferred.F;
Y_r = (A'*(Y- A*C - full(b)*f)) + C;
filtered.df = Y_r - repmat(residueF,1,size(C,2));
filtered.dfof = spdiags(filtered.F,0,K,K)\filtered.df;

% raw signal
CR = cell(K,2);    % store the information of the cell contour
rawA = A;          % spatial contour for raw signal
for idx=1:K
    A_temp = full(reshape(A(:,idx),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    [temp,ind] = sort(A_temp(:).^2,'ascend'); 
    temp =  cumsum(temp);
    ff = find(temp > (1-thr)*temp(end),1,'first');
    if ~isempty(ff)
        fp = find(A_temp >= A_temp(ind(ff)));
        [ii,jj] = ind2sub([d1,d2],fp);
        CR{idx,1} = [ii,jj]';
        CR{idx,2} = A_temp(fp)';
        rawA(fp,idx) = 1;
    end
end
nA = sqrt(sum(rawA.^2))';
rawA = rawA/spdiags(nA,0,K,K);    % normalize spatial components to unit energy

raw.df = rawA'*(Y-full(b)*f);     % no weighting on the ROI pixels, background substraction
residueF = getBaseline(raw.df,baselineRatio);
raw.F = median(rawA'*(full(b)*f),2) + residueF;
raw.df = raw.df - repmat(residueF,1,size(C,2));
raw.dfof = spdiags(raw.F,0,K,K)\raw.df;
raw.CR = CR;

end

% the following function get the baseline of an input signal
function [baseline] = getBaseline(signal,baselineRatio)
if nargin == 1
    baselineRatio = 0.25;
end
K = size(signal,1);
baseline = zeros(K,1);
for idx = 1:K
    temp=sort(signal(idx,:));
    index=find(temp>0);
    temp=temp(index);
    baseline(idx)=mean(temp(1:round(baselineRatio*length(temp))));
    if isnan(baseline(idx))
        baseline(idx) = 0;
    end
end
end