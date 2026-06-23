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


function pred_out = run_lsc_prediction_geoid_user(data_obs, pred_data, fit, ...
    distanceType, G, rho, muN, reg_rel)
% Run LSC prediction of geoid heights from gravity anomalies

if nargin < 8 || isempty(reg_rel)
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

pred = lsc_predict_geoid_from_gravity( ...
    lat_o, lon_o, g_o, std_o, ...
    lat_p, lon_p, ...
    fit, distanceType, G, rho, muN, reg_rel);

% ------ Output ------
pred_out.id      = id_p;
pred_out.lat     = lat_p;
pred_out.lon     = lon_p;
pred_out.h       = h_p;

pred_out.N_hat   = pred.N_hat;
pred_out.std_hat = pred.std_hat;
pred_out.var_hat = pred.var_hat;
pred_out.D_pred = pred.D_pred;

pred_out.fit = fit;
pred_out.distanceType = distanceType;
pred_out.G = G;
pred_out.rho = rho;
pred_out.muN = muN;
end