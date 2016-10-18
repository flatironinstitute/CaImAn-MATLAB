%% FOOPSI  (all done )
% AR1 
fprintf('AR1, FOOPSI...');
try
    ar1_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end

% AR2 
fprintf('\nAR2, FOOPSI...');
try
    ar2_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end

% kernel 
fprintf('\nkernel, FOOPSI...');
try
    kernel_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end


%% Constrained FOOPSI  (need AR2, kernel )
% AR1 
fprintf('\nAR1, constrained FOOPSI...');
try
    ar1_constrained_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end

%% Thresholded FOOPSI (all done)
% AR1 
fprintf('\nAR1, thresholded FOOPSI...');
try
    ar1_thresholded_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end

% AR2 
fprintf('\nAR2, thresholded FOOPSI... ');
try
    ar2_thresholded_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end


% AR2 
fprintf('\nkernel, thresholded FOOPSI... ');
try
    kernel_thresholded_foopsi;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end

%% MCMC 
% AR1 
fprintf('\nAR1, MCMC...');
try
    ar1_mcmc;
    fprintf('success!\n\n');
    drawnow; 
catch
    fprintf('fail!\n\n');
end
