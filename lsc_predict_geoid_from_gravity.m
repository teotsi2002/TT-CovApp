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


function pred = lsc_predict_geoid_from_gravity(lat_o, lon_o, g_o, std_o, ...
    lat_p, lon_p, fit, distanceType, G, rho, muN, reg_rel)
% LSC prediction of geoid height N from gravity anomalies using
% simplified cross-covariance model

if nargin < 12 || isempty(reg_rel)
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

if ~isscalar(G) || ~isscalar(rho) || ~isscalar(muN)
    error('G, rho and muN must be scalars.');
end

% Remove mean from gravity anomalies
mu_g = mean(g_o, 'omitnan');
gc = g_o - mu_g;

% Distances
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

% Auto-covariance of gravity anomalies
Cgg_oo = fit.model_fun(D_oo);
Cgg_po = fit.model_fun(D_po);
Cgg_pp = fit.model_fun(D_pp);

% Cross-covariance and auto-covariance for N
CNg_po = rho * G * Cgg_po;
CNN_pp = (G^2) * Cgg_pp;

% Observation noise
R = diag(std_o.^2);

% Regularized system matrix
lam = reg_rel * fit.sigma2;
A = Cgg_oo + R + lam * eye(No);
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

% Prediction of N
N_hat_centered = CNg_po * alpha;
N_hat = N_hat_centered + muN;

% Prediction error covariance
V = L \ CNg_po.';
Cov_pred = CNN_pp - V.' * V;

% Numerical symmetry enforcement
Cov_pred = 0.5 * (Cov_pred + Cov_pred.');

% Prediction error variances and standard deviations
var_hat = diag(Cov_pred);
var_hat(var_hat < 0) = 0;
std_hat = sqrt(var_hat);

% ------ Output ------
pred.N_hat   = N_hat;
pred.std_hat = std_hat;
pred.var_hat = var_hat;

pred.mu_g = mu_g;
pred.mu_N = muN;

pred.Coo = Cgg_oo;
pred.Cpo = CNg_po;
pred.Cpp = CNN_pp;
pred.D_pred = Cov_pred;

pred.rcondA = rcondA;
end