% ============================================================
% Author: Tsiridis Theodoros
% Year: 2026
%
% Diploma Thesis:
% Covariance Models in Point Adjustment Applications
% for the Analysis of Gravity Field Components
%
% Aristotle University of Thessaloniki (AUTH)
% Department of Surveying and Geodesy Engineering
%
% ------------------------------------------------------------
% Copyright (c) 2026 Tsiridis T.
% All rights reserved.
%
% This code is provided for academic purposes only.
% ============================================================


function [trainData, testData, idx_train, idx_test] = split_synthetic_data_for_lsc(data, train_fraction, seed)
%SPLIT_SYNTHETIC_DATA_FOR_LSC Split synthetic dataset into train/test sets

if nargin < 2 || isempty(train_fraction)
    train_fraction = 0.8;
end
if nargin < 3 || isempty(seed)
    seed = 1;
end

if train_fraction <= 0 || train_fraction >= 1
    error('train_fraction must be between 0 and 1.');
end

requiredFields = {'id','lat','lon','grav','std_grav'};
for k = 1:numel(requiredFields)
    if ~isfield(data, requiredFields{k})
        error('Input data must contain field "%s".', requiredFields{k});
    end
end

N = numel(data.id);
if N < 3
    error('Not enough points for train/test split.');
end

rng(seed);
idx = randperm(N);

Ntrain = round(train_fraction * N);
Ntrain = max(2, min(Ntrain, N-1));

idx_train = idx(1:Ntrain);
idx_test  = idx(Ntrain+1:end);

%  ------ Output ------
trainData.id       = data.id(idx_train);
trainData.lat      = data.lat(idx_train);
trainData.lon      = data.lon(idx_train);
trainData.grav     = data.grav(idx_train);
trainData.std_grav = data.std_grav(idx_train);

testData.id       = data.id(idx_test);
testData.lat      = data.lat(idx_test);
testData.lon      = data.lon(idx_test);
testData.grav     = data.grav(idx_test);
testData.std_grav = data.std_grav(idx_test);

if isfield(data, 'grav_true')
    trainData.grav_true = data.grav_true(idx_train);
    testData.grav_true  = data.grav_true(idx_test);
end
end