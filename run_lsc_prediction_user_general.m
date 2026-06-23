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


function pred_out = run_lsc_prediction_user_general(data_obs, pred_data, fit, distanceType, reg_rel)
%RUN_LSC_PREDICTION_USER_GENERAL Run LSC prediction for imported prediction points

if nargin < 5 || isempty(reg_rel)
    reg_rel = 1e-10;
end

lat_o = data_obs.lat(:);
lon_o = data_obs.lon(:);
g_o   = data_obs.grav(:);

if isfield(data_obs, 'std_grav') && ~isempty(data_obs.std_grav)
    std_o = data_obs.std_grav(:);
else
    error('Observation dataset must include std_grav for LSC.');
end

id_p  = pred_data.id(:);
lat_p = pred_data.lat(:);
lon_p = pred_data.lon(:);

if isfield(pred_data, 'h') && ~isempty(pred_data.h)
    h_p = pred_data.h(:);
else
    h_p = nan(size(lat_p));
end

lat_min = min(lat_o); lat_max = max(lat_o);
lon_min = min(lon_o); lon_max = max(lon_o);

if any(lat_p < lat_min | lat_p > lat_max | lon_p < lon_min | lon_p > lon_max)
    warning(['Some prediction points lie outside the observation bounding box. ', ...
        'This may reduce prediction reliability.']);
end

pred = lsc_predict_general(lat_o, lon_o, g_o, std_o, lat_p, lon_p, fit, distanceType, reg_rel);

% ------ Output ------
pred_out.id      = id_p;
pred_out.lat     = lat_p;
pred_out.lon     = lon_p;
pred_out.h       = h_p;

pred_out.g_hat   = pred.g_hat;
pred_out.std_hat = pred.std_hat;
pred_out.var_hat = pred.var_hat;
pred_out.D_pred = pred.D_pred;
pred_out.fit = fit;
pred_out.distanceType = distanceType;
end