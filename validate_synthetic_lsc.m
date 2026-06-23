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


function val = validate_synthetic_lsc(data, modelID, train_fraction, bin_width, max_dist, reg_rel)
% VALIDATE_SYNTHETIC_LSC Validate LSC on synthetic field with train/test split

if nargin < 3 || isempty(train_fraction)
    train_fraction = 0.8;
end
if nargin < 4
    bin_width = [];
end
if nargin < 5    max_dist = [];
end
if nargin < 6 || isempty(reg_rel)
    reg_rel = 1e-10;
end

% ------ Basic checks ------
if ~isfield(data, 'lat') || ~isfield(data, 'lon') || ~isfield(data, 'grav') || ~isfield(data, 'std_grav')
    error('data must contain lat, lon, grav, std_grav.');
end

N = numel(data.grav);
if N < 10
    error('Not enough points for validation.');
end

% ------ Random split ------
rng(1);  % reproducible split
idx = randperm(N);

Ntrain = round(train_fraction * N);
Ntrain = max(2, min(Ntrain, N-1));

idx_train = idx(1:Ntrain);
idx_test  = idx(Ntrain+1:end);

% ------ Training set ------
lat_o  = data.lat(idx_train);
lon_o  = data.lon(idx_train);
g_o    = data.grav(idx_train);
std_o  = data.std_grav(idx_train);

% ------ Prediction/test set ------
lat_p = data.lat(idx_test);
lon_p = data.lon(idx_test);

g_test_obs = data.grav(idx_test);

has_true = isfield(data, 'grav_true') && ~isempty(data.grav_true);
if has_true
    g_test_true = data.grav_true(idx_test);
else
    g_test_true = [];
end

% ------ Empirical covariance from training set ------
results_emp = Empirical_Model(lat_o, lon_o, g_o, bin_width, max_dist);

% ------ Dynamic model fitting ------
fit = dynamic_model_fitting_test(results_emp, modelID);

% ------ LSC prediction ------
pred = lsc_predict_auto(lat_o, lon_o, g_o, std_o, lat_p, lon_p, fit, reg_rel);

% ------ Error wrt observed test values ------
err_obs = pred.g_hat - g_test_obs;
rmse_obs = sqrt(mean(err_obs.^2, 'omitnan'));
mae_obs  = mean(abs(err_obs), 'omitnan');
bias_obs = mean(err_obs, 'omitnan');

% ------ Error wrt true synthetic field ------
if has_true
    err_true = pred.g_hat - g_test_true;
    rmse_true = sqrt(mean(err_true.^2, 'omitnan'));
    mae_true  = mean(abs(err_true), 'omitnan');
    bias_true = mean(err_true, 'omitnan');
else
    err_true = [];
    rmse_true = NaN;
    mae_true  = NaN;
    bias_true = NaN;
end

% ------ Output ------
val.idx_train = idx_train;
val.idx_test  = idx_test;

val.results_emp = results_emp;
val.fit  = fit;
val.pred = pred;

val.err_obs  = err_obs;
val.rmse_obs = rmse_obs;
val.mae_obs  = mae_obs;
val.bias_obs = bias_obs;

val.err_true  = err_true;
val.rmse_true = rmse_true;
val.mae_true  = mae_true;
val.bias_true = bias_true;
end