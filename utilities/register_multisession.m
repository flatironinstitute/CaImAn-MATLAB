function [A_union, assignments, matchings] = register_multisession(A, options, templates, options_mc)
% REGISTER_MULTISESSION - register ROIs from multiple recording sessions
%
%   [A_UNION, ASSIGNMENTS, MATCHINGS] = REGISTER_MULTISESSION(A,...
%          OPTIONS, TEMPLATES, OPTIONS_MC)
%
% Register ROIs across multiple sessions using an intersection over union metric
% and the Hungarian algorithm for optimal matching. Registration occurs by 
% aligning session 1 to session 2, keeping the union of the matched and 
% non-matched components to register with session 3 and so on.
%
% INPUTS:
% A:                      cell array with spatial components from each session
% options                 parameter structure with inputs:
%       d1:               number of rows in FOV
%       d2:               number of columns in FOV
%       d3:               number of planes in FOV (default: 1)
%       dist_maxthr:      threshold for turning spatial components into binary masks (default: 0.1)
%       dist_exp:         power n for distance between masked components: dist = 1 - (and(m1,m2)/or(m1,m2))^n (default: 1)
%       dist_thr:         threshold for setting a distance to infinity. (default: 0.5)
%       dist_overlap_thr: overlap threshold for detecting if one ROI is a subset of another (default: 0.8)
%       template1:        template from motion correction of the first session
%       template2:        template from motion correction of the second session
%       options_mc:       motion correction options
%       plot_reg:         create a contour plot of registered ROIs
% templates:              cell array with a reference image from each session
% options_mc:             normcorre options structure (for template alignment) 

% OUTPUTS:
% A_union:                union of ROIs aligned to session 1 (for matched
%                               pairs the ROIs from session # 1 are kept)
% assignments:            matrix of size # of total distinct components x # sessions
%                               element [i,j] = k if component k from session j is mapped to component
%                               i in the A_union matrix. If there is no much the value is NaN
% matchings:              cell array. Each entry matchings{i}(j) = k means that component j from session
%                               i is represented by component k in A_union

defoptions = CNMFSetParms;
if ~exist('options','var'); options = defoptions; end

if ~exist('templates','var') || ~exist('options_mc','var');
    warning('Some required inputs for aligning ROIs before registering are missing. Skipping alignment');
    align_flag = false;
else
    align_flag = true;
end

if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); options.d1 = d1; end 
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); options.d2 = d2; end 
if ~isfield(options,'d3') || isempty(options.d3); options.d3 = 1; end
if ~isfield(options,'dist_maxthr') || isempty(options.maxthr); options.dist_maxthr = 0.15; end
if ~isfield(options,'dist_exp') || isempty(options.dist_exp); options.dist_exp = 1; end
if ~isfield(options,'dist_thr') || isempty(options.dist_thr); options.dist_thr = 0.5; end
if ~isfield(options,'dist_overlap_thr') || isempty(options.dist_overlap_thr); options.dist_overlap_thr = 0.8; end
options.plot_reg = false;

n_sessions = length(A);
if ~align_flag
    templates = {[]};
end
if length(templates) == 1
    for i = 2:n_sessions
        templates{i} = templates{1};
    end
end

siz = [options.d1,options.d2,options.d3];
options_mc.correct_bidir = false;
A_union = A{1};
matchings{1} = 1:size(A{1},2);

for session = 2:n_sessions
   [matched_ROIs,nonmatched_1,nonmatched_2,A2,R,A_un] = register_ROIs(A{session},A_union,options,templates{session},templates{session-1},options_mc);
   A_union = A_un;
   A_union(:, matched_ROIs(:,2)) = A{session}(:, matched_ROIs(:,1));
   A_union = [A_union, A{session}(:,nonmatched_1)];
   new_match = zeros(1,size(A{session},2));
   new_match(matched_ROIs(:,1)) = matched_ROIs(:,2);
   new_match(nonmatched_1) = size(A_un,2)+1:size(A_union,2);
   matchings{session} = new_match;
end

assignments = NaN(size(A_union,2), n_sessions);
for session = 1:n_sessions
    assignments(matchings{session}, session) = 1:length(matchings{session});
end