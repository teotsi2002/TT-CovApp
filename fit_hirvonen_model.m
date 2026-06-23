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


function fit = fit_hirvonen_model(results)
%FIT_HIRVONEN_MODEL Fit Hirvonen covariance model to empirical covariance

% Use zero-lag + valid empirical bins
h     = [0; results.lags_valid(:)];
C_emp = [results.C0_emp; results.C_emp_valid(:)];
w     = [max(results.pair_count_valid); results.pair_count_valid(:)];

if isempty(h) || isempty(C_emp)
    error('No valid empirical covariance points available for fitting.');
end

% Initial guesses
sigma2_0 = max(results.C0_emp, eps);

if isfinite(results.Lhalf_emp) && results.Lhalf_emp > 0
    a_0 = results.Lhalf_emp;
else
    a_0 = max(h) / 3;
end

a_0 = max(a_0, eps);

% Work in log-space
theta0 = [log(sigma2_0); log(a_0)];

% Weighted SSE objective
obj = @(th) sum( w .* ...
    (C_emp - cov_hirvonen(h, exp(th(1)), exp(th(2)))).^2 );

% Minimize
options = optimset('Display','off', ...
    'TolX',1e-8, ...
    'TolFun',1e-8, ...
    'MaxIter',1000, ...
    'MaxFunEvals',4000);

theta_hat = fminsearch(obj, theta0, options);

% Estimated parameters
sigma2_hat = exp(theta_hat(1));
a_hat      = exp(theta_hat(2));

% Fitted covariance on empirical lags
C_fit = cov_hirvonen(h, sigma2_hat, a_hat);

% Final SSE
sse = obj(theta_hat);

% Output
fit.sigma2    = sigma2_hat;
fit.a         = a_hat;
fit.theta_hat = theta_hat;
fit.C_fit     = C_fit;
fit.sse0      = obj(theta0);
fit.sse       = sse;
fit.Lhalf_mod = a_hat;

fit.model_name = 'Hirvonen';
fit.model_fun  = @(h) cov_hirvonen(h, sigma2_hat, a_hat);
end

% COV_HIRVONEN Hirvonen covariance model
function C = cov_hirvonen(h, sigma2, a)

C = sigma2 ./ (1 + (h./a).^2);
end