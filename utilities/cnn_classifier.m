function [ind,value] = cnn_classifier(A,dims,classifier,thr)

%cnn_classifer classify spatial components using a pretrained CNN
%classifier using the keras importer add on.
%   IND = cnn_classifier(A,dims,classifier,thr) returns a binary vector indicating
%   whether the set of spatial components A, with dimensions of the field
%   of view DIMS, pass the threshold THR for the given CLASSIFIER
%
%   [IND,VALUE] = cnn_classifier(A,dims,classifier,thr) also returns the
%   output value of the classifier 
%
%   INPUTS:
%   A:              2d matrix
%   dims:           vector with dimensions of the FOV
%   classifier:     path to pretrained classifier model (downloaded if it
%                       doesn't exist)
%   thr:            threshold for accepting component (default: 0.2)
%
%   note: The function requires Matlab version 2017b (9.3) or later, Neural
%   Networks toolbox version 2017b (11.0) or later, the Neural Network 
%   Toolbox(TM) Importer for TensorFlow-Keras Models.

%   Written by Eftychios A. Pnevmatikakis. Classifier trained by Andrea
%   Giovannucci, Flatiron Institute, 2017

K = size(A,2);                          % number of components

if verLessThan('matlab','9.3') || verLessThan('nnet','11.0') || isempty(which('importKerasNetwork'))
    warning(strcat('The function cnn_classifier requires Matlab version 2017b (9.3) or later, Neural\n', ...
        'Networks toolbox version 2017b (11.0) or later, the Neural Networks ', ...
        'Toolbox(TM) Importer for TensorFlow-Keras Models.'))
    ind = true(K,1);
    value = ones(K,1);
else
    if ~exist('thr','var'); thr = 0.2; end

    A = A/spdiags(sqrt(sum(A.^2,1))'+eps,0,K,K);      % normalize to sum 1 for each compoennt
    A_com = extract_patch(A,dims,[50,50]);  % extract 50 x 50 patches

    if ~exist(classifier,'file')
        url = 'https://www.dropbox.com/s/1csymylbne7yyt0/cnn_model.h5?dl=1';
        classifier = 'cnn_model.h5';
        outfilename = websave(classifier,url);
    end

    net_classifier = importKerasNetwork(classifier,'ClassNames',["rejected","accepted"]);
    out = predict(net_classifier,double(A_com));
    value = out(:,2);
    ind = (value >= thr);
end