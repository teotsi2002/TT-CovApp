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


function pred = lsc_predict_general(lat_o, lon_o, g_o, std_o, lat_p, lon_p, fit, distanceType, reg_rel)
% LSC_PREDICT_GENERAL Least Squares Collocation prediction with selectable horizontal distances

if nargin < 9 || isempty(reg_rel)
    reg_rel = 1e-10;
end

lat_o = lat_o(:);
lon_o = lon_o(:);
g_o   = g_o(:);
std_o = std_o(:);

lat_p = lat_p(:);
lon_p = lon_p(:);

No = numel(g_o);
Np = numel(lat_p);

if numel(lat_o) ~= No || numel(lon_o) ~= No || numel(std_o) ~= No
    error('Observation vectors lat_o, lon_o, g_o, std_o must have the same length.');
end

if numel(lon_p) ~= Np
    error('Prediction coordinate vectors lat_p and lon_p must have the same length.');
end

if isempty(fit) || ~isfield(fit, 'model_fun') || ~isfield(fit, 'sigma2')
    error('fit structure must contain fit.model_fun and fit.sigma2.');
end

mu = mean(g_o, 'omitnan');
gc = g_o - mu;

switch distanceType
    case 'Spherical'
        D_oo = compute_sph_distance_LSC(lat_o, lon_o, lat_o, lon_o);
        D_po = compute_sph_distance_LSC(lat_p, lon_p, lat_o, lon_o);
        D_pp = compute_sph_distance_LSC(lat_p, lon_p, lat_p, lon_p);

    case 'Euclidean'
        D_oo = compute_planar_distance_LSC(lat_o, lon_o, lat_o, lon_o);
        D_po = compute_planar_distance_LSC(lat_p, lon_p, lat_o, lon_o);
        D_pp = compute_planar_distance_LSC(lat_p, lon_p, lat_p, lon_p);

    otherwise
        error('Unknown distanceType. Use "Spherical" or "Euclidean".');
end

Coo = fit.model_fun(D_oo);
Cpo = fit.model_fun(D_po);
Cpp = fit.model_fun(D_pp);

R = diag(std_o.^2);
lam = reg_rel * fit.sigma2;

A = Coo + R + lam * eye(No);
A = (A + A.') / 2;

rcondA = rcond(A);

[L, p] = chol(A, 'lower');
if p > 0
    lam_boost = 1e-6 * fit.sigma2;
    A = A + lam_boost * eye(No);
    [L, p] = chol(A, 'lower');
    if p > 0
        error('Cholesky failed even after regularization boost.');
    end
end

alpha = L' \ (L \ gc);

g_hat_centered = Cpo * alpha;
g_hat = g_hat_centered + mu;

V = L \ Cpo.';
Cov_pred = Cpp - V.' * V;

% Numerical symmetry enforcement
Cov_pred = 0.5 * (Cov_pred + Cov_pred.');

var_hat = diag(Cov_pred);
var_hat(var_hat < 0) = 0;
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