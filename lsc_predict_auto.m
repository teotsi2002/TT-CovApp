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


function pred = lsc_predict_auto(lat_o, lon_o, g_o, std_o, lat_p, lon_p, fit, reg_rel)
%LSC_PREDICT_AUTO Least Squares Collocation prediction using auto-covariance

if nargin < 8 || isempty(reg_rel)
    reg_rel = 1e-10;
end

% ------ Ensure column vectors ------
lat_o = lat_o(:);
lon_o = lon_o(:);
g_o   = g_o(:);
std_o = std_o(:);

lat_p = lat_p(:);
lon_p = lon_p(:);

No = numel(g_o);
Np = numel(lat_p);

% ------ Input checks ------
if numel(lat_o) ~= No || numel(lon_o) ~= No || numel(std_o) ~= No
    error('Observation vectors lat_o, lon_o, g_o, std_o must have the same length.');
end

if numel(lon_p) ~= Np
    error('Prediction coordinate vectors lat_p and lon_p must have the same length.');
end

if isempty(fit) || ~isfield(fit, 'model_fun') || ~isfield(fit, 'sigma2')
    error('fit structure must contain at least fit.model_fun and fit.sigma2.');
end

% ------ Remove mean from observations ------
mu = mean(g_o, 'omitnan');
gc = g_o - mu;

% ------ Distance matrices ------
D_oo = compute_sph_distance_LSC(lat_o, lon_o, lat_o, lon_o);
D_po = compute_sph_distance_LSC(lat_p, lon_p, lat_o, lon_o);
D_pp = compute_sph_distance_LSC(lat_p, lon_p, lat_p, lon_p);

% ------ Covariance matrices from fitted model ------
Coo = fit.model_fun(D_oo);
Cpo = fit.model_fun(D_po);
Cpp = fit.model_fun(D_pp);

% ------ Observation noise covariance ------
R = diag(std_o.^2);

% ------ Regularization ------
lam = reg_rel * fit.sigma2;

% ------ System matrix ------
A = Coo + R + lam * eye(No);
A = (A + A.') / 2;   % enforce symmetry

% ------ Conditioning check ------
rcondA = rcond(A);

% ------ Cholesky factorization with fallback ------
[L, p] = chol(A, 'lower');

if p > 0
    lam_boost = 1e-6 * fit.sigma2;
    A = A + lam_boost * eye(No);
    [L, p] = chol(A, 'lower');

    if p > 0
        error('Cholesky failed even after regularization boost.');
    end
end

% ------ Solve for collocation weights ------
alpha = L' \ (L \ gc);

% ------ Prediction ------
g_hat_centered = Cpo * alpha;
g_hat = g_hat_centered + mu;

% ------ Prediction error covariance ------
V = L \ Cpo.';
Cov_pred = Cpp - V.' * V;

% Numerical symmetry enforcement
Cov_pred = 0.5 * (Cov_pred + Cov_pred.');

% Prediction error variances and standard deviations
var_hat = diag(Cov_pred);
var_hat(var_hat < 0) = 0;  % numerical safety
std_hat = sqrt(var_hat);

% ------ Output ------
pred.g_hat   = g_hat;
pred.std_hat = std_hat;
pred.var_hat = var_hat;
pred.mu      = mu;

pred.Coo = Coo;
pred.Cpo = Cpo;
pred.Cpp = Cpp;
pred.D_pred = Cov_pred;

pred.rcondA = rcondA;
end