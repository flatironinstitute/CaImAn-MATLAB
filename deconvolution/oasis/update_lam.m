function [c, active_set, lam, s, flag_lam] = update_lam(y, c, active_set, g, lam, thresh)
%% update the tuning parameter lambda to reduce |s|_1 while |y-c|_2^2 <= sn^2*T

%% inputs:
%   y:  T*1 vector, One dimensional array containing the fluorescence intensities 
        %withone entry per time-bin.
%   c: T*1 vector, previous c 
%   active_set: npools*4 matrix, previous active sets 
%   g:  scalar, Parameter of the AR(1) process that models the fluorescence ...
        %impulse response.   
%   lam:  scalar, curret value of sparsity penalty parameter lambda. 
%   thresh: maximum residual sn*T^2
       
%% outputs
%   c: T*1 vector
%   s: T*1 vector, spike train 
%   active_set: npool x 4 matrix, active sets 
%   lam: scalar, new tuning parameter 
%   flag_lam: bool, True if it update lam with success
%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% ported from the Python implementation from Johannes Friedrich

%% References 
% Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging

c = reshape(c, [],1); 
y = reshape(y, [],1); 
len_active_set = size(active_set,1); 

res = y - c; 
RSS = res'*res; 

temp = zeros(size(c)); 
for ii=1:len_active_set
    ti = active_set(ii, 3); 
    li = active_set(ii, 4); 
    idx = 0:(li-1); 
    temp(ti+idx) = (1-g^li)/ active_set(ii,2) * g.^(idx);   
end

aa = temp'*temp; 
bb = res'*temp; 
cc = RSS-thresh; 
ll = (-bb + sqrt(bb^2-aa*cc)) / aa; 
if imag(ll)~=0
    flag_lam = false; 
    s = [0; c(2:end)-c(1:(end-1))*g]; 
    return; 
else
    flag_lam = true; 
end
lam = lam + ll; 

active_set(:,1) = active_set(:,1) - ll*(1-g.^active_set(:,4)); 
[c, s, active_set] = oasisAR1(y, g, lam, 0, active_set); 























