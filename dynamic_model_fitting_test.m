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


function fit = dynamic_model_fitting_test(results, modelID)
% FIT_COVARIANCE_MODEL Dynamic fitting of covariance model

% Special case: Forsberg planar (3-parameter model)

if strcmpi(modelID, 'forsbergplanar')
    fit = fit_forsberg_planar_model(results);
    return
end


% Common fitting logic for 2-parameter models


% Use zero-lag + valid bins
h     = [0; results.lags_valid(:)];
C_emp = [results.C0_emp; results.C_emp_valid(:)];
w     = [max(results.pair_count_valid); results.pair_count_valid(:)];

if isempty(h) || isempty(C_emp)
    error('No valid empirical covariance points available for fitting.');
end

sigma2_0 = max(results.C0_emp, eps);

switch lower(modelID)

    case 'exponential'
        if isfinite(results.Lhalf_emp) && results.Lhalf_emp > 0
            a_0 = results.Lhalf_emp / log(2);
        else
            a_0 = max(h) / 3;
        end

        model_fun_raw = @(hh, s2, a) s2 .* exp(-hh ./ a);
        Lhalf_fun = @(a) a * log(2);
        model_name = 'Exponential';

    case 'gaussian'
        if isfinite(results.Lhalf_emp) && results.Lhalf_emp > 0
            a_0 = results.Lhalf_emp / sqrt(log(2));
        else
            a_0 = max(h) / 3;
        end

        model_fun_raw = @(hh, s2, a) s2 .* exp(-(hh.^2) ./ (a.^2));
        Lhalf_fun = @(a) a * sqrt(log(2));
        model_name = 'Gaussian';

    case 'hirvonen'
        if isfinite(results.Lhalf_emp) && results.Lhalf_emp > 0
            a_0 = results.Lhalf_emp;
        else
            a_0 = max(h) / 3;
        end

        model_fun_raw = @(hh, s2, a) s2 ./ (1 + (hh ./ a).^2);
        Lhalf_fun = @(a) a;
        model_name = 'Hirvonen';

    otherwise
        error('Unknown modelID: %s', modelID);
end

a_0 = max(a_0, eps);

% Work in log-space
theta0 = [log(sigma2_0); log(a_0)];

% Weighted SSE objective
obj = @(th) sum(w .* (C_emp - model_fun_raw(h, exp(th(1)), exp(th(2)))).^2);

% Optimizer options
opts = optimset('Display','off', ...
    'TolX',1e-8, ...
    'TolFun',1e-8, ...
    'MaxIter',1000, ...
    'MaxFunEvals',4000);

theta_hat = fminsearch(obj, theta0, opts);

sigma2_hat = exp(theta_hat(1));
a_hat      = exp(theta_hat(2));

C_fit = model_fun_raw(h, sigma2_hat, a_hat);

% ----- Output ------
fit.modelID    = lower(modelID);
fit.model_name = model_name;

fit.sigma2     = sigma2_hat;
fit.a          = a_hat;
fit.theta_hat  = theta_hat;

fit.C_fit      = C_fit;
fit.sse0       = obj(theta0);
fit.sse        = obj(theta_hat);

fit.Lhalf_mod  = Lhalf_fun(a_hat);
fit.model_fun  = @(hh) model_fun_raw(hh, sigma2_hat, a_hat);
end