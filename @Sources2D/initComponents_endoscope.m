function [center, Cn, PNR] = initComponents_endoscope(obj, Y, K, patch_sz, debug_on, save_avi)
%% initializing spatial/temporal components for miceoendoscopic data
%% input:
%   Y:  d1 X d2 X T matrix or (d1*d2)X T matrix, video data
%   K:  scalar, maximum number of neurons
%   patch_sz: scalar or vector, if its scalar, then the patch size is
%       patch_sz * patch_sz; if its vector, then the optical field is divided
%       into patch_sz(1) X patch_sz(2) patches.
%   debug_on: options for showing procedure of detecting neurons
%   save_avi: save the video of initialization procedure
%% Output:
%   center: d*2 matrix, centers of all initialized neurons.
%   Cn:     correlation image
%   PNR:    peak to noise ratio
%% Author: Pengcheng Zhou, Carnegie Mellon University, zhoupc1988@gmail.com

%% process parameters
d1 = obj.options.d1;
d2 = obj.options.d2;

if isfield(obj.options, 'nk') % number of knots for creating spline basis
    nk = obj.options.nk;
else
    nk = 1;
end
% maximum neuron number in each patch
if (~exist('K', 'var')) || (isempty(K))
    % if K is not specified, use a very large number as default
    K = round((d1*d2));
end

% exporting initialization procedures as a video
if ~exist('save_avi', 'var')||isempty(save_avi)
    save_avi = false;
elseif ischar(save_avi) || (save_avi==true)
    debug_on = true; % turn on debug mode
else
    save_avi = false; %don't save initialization procedure
end

% debug mode and exporting results
if ~exist('debug_on', 'var')
    debug_on = false;
end

options = obj.options;
options.kernel = obj.kernel;
% divide the optical fields into multiple patches
if (~exist('patch_sz', 'var'))||(isempty(patch_sz))||(max(patch_sz(:))==1)
    % use the whole optical field directly
    if nk>1  % detrend data
        Ydt = detrend_data(obj.reshape(double(Y), 1), nk); % detrend data
        [results, center, Cn, PNR] = greedyROI_endoscope(Ydt, K, options, debug_on, save_avi);
    else    % without detrending
        [results, center, Cn, PNR] = greedyROI_endoscope(Y, K, options, debug_on, save_avi);
    end
    obj.A = results.Ain;
    obj.C = results.Cin;
    obj.C_raw = results.Cin_raw;
    if options.deconv_flag
        obj.S = results.Sin;
        obj.P.kernel_pars = results.kernel_pars;
    else
        obj.S = zeros(size(obj.C));
    end
    obj.Cn = Cn;
    return;
elseif isscalar(patch_sz)
    % patch size has been assigned
    r0_patch = (1:patch_sz:d1);     % break points in the y axis
    r0_patch(end) = d1;             % break points in the x axis
    c0_patch = (1:patch_sz:d2);
    c0_patch(end) = d2;
elseif length(patch_sz)==2
    % the patch number has been specified.
    r0_patch = round(linspace(1, d1, 1+patch_sz(1)));
    c0_patch = round(linspace(1, d2, 1+patch_sz(2)));
else
    fprintf('wrong parameters for patch_sz\n\n');
    return;
end

%% start initialization
nr_patch = length(r0_patch)-1;  % number of patches in y axis
nc_patch = length(c0_patch)-1;  % number of patches in x axis
Y = obj.reshape(Y, 2);   % represent each frame with a matrix
if ~isfield(options, 'bd') || isempty(options.bd')
    options.bd = options.gSiz;   % boundary pixesl to be ignored during the process of detecting seed pixels
end
bd = options.bd;
Ain = cell(nr_patch, nc_patch); % save spatial components of neurons in each patch
Cin = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, denoised trace
Sin = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, deconvolved trace
Cin_raw = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch, raw trace
kernel_pars = cell(nr_patch, nc_patch); % save temporal components of neurons in each patch
center = cell(nr_patch, nc_patch);     % save centers of all initialized neurons
Cn = zeros(d1, d2);
PNR = zeros(d1, d2);

% initialize A and C for each patch
flag_patch = false(nr_patch, nc_patch);
for mr = 1:nr_patch
    r0 = max(1, r0_patch(mr)-bd); % minimum row index of the patch
    r1 = min(d1, r0_patch(mr+1)+bd-1); % maximum row index of the patch
    for mc = 1:nc_patch
        c0 = max(1, c0_patch(mc)-bd); % minimum column index of the patch
        c1 = min(d2, c0_patch(mc+1)+bd-1); % maximum column index of the patch
        tmp_options = options;
        tmp_options.d1 = (r1-r0)+1;
        tmp_options.d2 = (c1-c0)+1;
        tmp_options.kernel = obj.kernel;
        
        % name of the video
        if ischar(save_avi)
            tmp_save_avi = sprintf('%s_%d_%d_%d_%d.avi', save_avi, r0, c0, r1, c1);
        else
            tmp_save_avi = save_avi;
        end
        
        % take the patch from the raw data
        nrows = (r1-r0+1);  % number of rows in the patch
        ncols = (c1-c0+1);  %number of columns in the patch
        Ypatch = double(reshape(Y(r0:r1, c0:c1, :), nrows*ncols, []));
        
        % top patch, some signals have been explained by neurons in the top
        % patch
        if mr>1
            tmpA = reshape(Ain{mr-1, mc}, d1, d2, []);
            tmpC = Cin{mr-1, mc};
            Ypatch = Ypatch - reshape(tmpA(r0:r1, c0:c1,:), nrows*ncols, [])*tmpC;
        end
        
        % left patch, some signals have been explained by neurons in the
        % left patch
        if mc>1
            tmpA = reshape(Ain{mr, mc-1}, d1, d2, []);
            tmpC = Cin{mr, mc-1};
            Ypatch = Ypatch-reshape(tmpA(r0:r1, c0:c1,:), nrows*ncols, [])*tmpC;
        end
        if nk>1
            Ypatch_dt = detrend_data(Ypatch, nk); % detrend data
            [tmp_results, tmp_center, tmp_Cn, tmp_PNR, save_avi] = greedyROI_endoscope(Ypatch_dt, K, tmp_options, debug_on, tmp_save_avi);
        else
            [tmp_results, tmp_center, tmp_Cn, tmp_PNR, save_avi] = greedyROI_endoscope(Ypatch, K, tmp_options, debug_on, tmp_save_avi);
        end
        close(gcf);
        tmp_Ain = tmp_results.Ain;
        tmp_Cin = tmp_results.Cin;
        tmp_Sin = tmp_results.Sin;
        tmp_Cin_raw = tmp_results.Cin_raw;
        tmp_kernel_pars = tmp_results.kernel_pars;
        tmp_K = size(tmp_Ain, 2);   % number of neurons within the selected patch
        temp = zeros(d1, d2, tmp_K);  % spatial components of all neurosn
        temp(r0:r1, c0:c1, :) = reshape(full(tmp_Ain), tmp_options.d1, tmp_options.d2, []);
        Ain{mr, mc} = reshape(temp, d1*d2, tmp_K);
        Cin{mr, mc} = tmp_Cin;      % temporal components of all neurons
        Cin_raw{mr, mc} = tmp_Cin_raw;
        if tmp_options.deconv_flag
            Sin{mr, mc} = tmp_Sin;
            kernel_pars{mr,mc} = tmp_kernel_pars;
        end
        center{mr, mc} = bsxfun(@plus, tmp_center, [r0-1, c0-1]); % centers
        Cn(r0:r1, c0:c1) = max(Cn(r0:r1, c0:c1), tmp_Cn);
        PNR(r0:r1, c0:c1) = max(PNR(r0:r1, c0:c1), tmp_PNR);
        
        % display initialization progress
        flag_patch(mr, mc) = true;
        clc;
        fprintf('\n Progress......\n');
        for tmpm = 1:nr_patch
            for tmpn=1:nc_patch
                if flag_patch(tmpm, tmpn)
                    fprintf('X\t');
                else
                    fprintf('O\t');
                end
            end
            fprintf('\n');
        end
    end
end

%% export the results
Ain = cell2mat(reshape(Ain, 1, []));
Cin = cell2mat(reshape(Cin, [], 1));
Cin_raw = cell2mat(reshape(Cin_raw, [], 1));
center = cell2mat(reshape(center, [], 1));
obj.A = Ain;
obj.C = Cin;
obj.C_raw = Cin_raw;
if tmp_options.deconv_flag
    Sin = cell2mat(reshape(Sin, [], 1));
    kernel_pars = cell2mat(reshape(kernel_pars, [], 1));
    obj.S = Sin;
    obj.P.kernel_pars = kernel_pars;
else
    obj.S = zeros(size(obj.C));
end
obj.Cn = Cn;
end