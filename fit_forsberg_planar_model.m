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


function fit = fit_forsberg_planar_model(results)
% FIT_FORSBERG_PLANAR_MODEL Fit Forsberg planar covariance model

if ~isfield(results, 'distance_type') || ~strcmpi(results.distance_type, 'planar')
    error('Forsberg planar fitting requires planar empirical lags.');
end

% Use zero-lag + valid bins
h     = [0; results.lags_valid(:)];
C_emp = [results.C0_emp; results.C_emp_valid(:)];
w     = [max(results.pair_count_valid); results.pair_count_valid(:)];

if isempty(h) || isempty(C_emp)
    error('No valid empirical covariance points available for fitting.');
end

% Initial guesses
C0_0 = max(results.C0_emp, eps);

if isfinite(results.Lhalf_emp) && results.Lhalf_emp > 0
    D_0 = max(results.Lhalf_emp / 4, eps);
else
    D_0 = max(max(h) / 10, eps);
end

T_0 = max(3 * D_0, eps);

% Parameterization: theta = [log(C0); log(D); log(T)]
theta0 = [log(C0_0); log(D_0); log(T_0)];

% Weighted SSE objective
obj = @(th) sum(w .* ...
    (C_emp - cov_forsberg_planar(h, 0, 0, exp(th(1)), exp(th(2)), exp(th(3)))).^2);

% Optimizer options
opts = optimset('Display','off', ...
    'TolX',1e-6, ...
    'TolFun',1e-6, ...
    'MaxIter',800, ...
    'MaxFunEvals',3000);

theta_hat = fminsearch(obj, theta0, opts);

% Estimated parameters
C0_hat = exp(theta_hat(1));
D_hat  = exp(theta_hat(2));
T_hat  = exp(theta_hat(3));

% Fitted covariance on empirical lags
C_fit = cov_forsberg_planar(h, 0, 0, C0_hat, D_hat, T_hat);

% Final SSE
sse = obj(theta_hat);

% Numerical half-covariance distance
h_plot = linspace(0, max(h), 1000).';
C_plot = cov_forsberg_planar(h_plot, 0, 0, C0_hat, D_hat, T_hat);

target = 0.5 * C0_hat;
Lhalf_mod = NaN;

k = find(C_plot <= target, 1, 'first');
if ~isempty(k) && k > 1
    x1 = h_plot(k-1); x2 = h_plot(k);
    y1 = C_plot(k-1); y2 = C_plot(k);

    if y2 ~= y1
        Lhalf_mod = x1 + (target - y1) * (x2 - x1) / (y2 - y1);
    end
end

% ------ Output ------
fit.modelID    = 'forsbergplanar';
fit.model_name = 'Forsberg planar';

fit.C0         = C0_hat;
fit.D          = D_hat;
fit.T          = T_hat;
fit.theta_hat  = theta_hat;

fit.sigma2     = C0_hat;
fit.a          = D_hat;

fit.C_fit      = C_fit;
fit.sse0       = obj(theta0);
fit.sse        = sse;
fit.Lhalf_mod  = Lhalf_mod;

fit.model_fun  = @(hh) cov_forsberg_planar(hh, 0, 0, C0_hat, D_hat, T_hat);
end