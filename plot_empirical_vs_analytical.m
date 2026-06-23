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


function plot_empirical_vs_analytical(results, fit)
%PLOT_EMPIRICAL_VS_MODEL Plot empirical covariance and fitted analytical model.
%
% INPUT:
%   results : structure returned by Empirical_Model
%   fit     : structure returned by fitting function (e.g. fit_exponential_model)

% Extract empirical data
lags_all   = results.lags;
C_all      = results.C_emp;
lags_valid = results.lags_valid;
C_valid    = results.C_emp_valid;
C0_emp     = results.C0_emp;
Lhalf_emp  = results.Lhalf_emp;
min_pairs  = results.min_pairs_bin;

% Dense grid for smooth model curve
h_plot = linspace(0, max(lags_valid, [], 'omitnan'), 500).';
C_plot = fit.model_fun(h_plot);

figure('Color','w');
hold on;
grid on;
box on;

% Empirical covariance (all bins)
h1 = plot(lags_all, C_all, '-o', ...
    'LineWidth', 1.0, ...
    'MarkerSize', 4);

% Valid bins
h2 = plot(lags_valid, C_valid, 'o', ...
    'LineWidth', 1.5, ...
    'MarkerSize', 5);

% Analytical fitted model
h3 = plot(h_plot, C_plot, 'LineWidth', 2.0);

% Half-covariance level for empirical
target_emp = 0.5 * C0_emp;
h4 = yline(target_emp, '--', '$0.5\,C(0)$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','left');

% Optional: model half-covariance line
if isfield(fit, 'sigma2')
    target_mod = 0.5 * fit.sigma2;
    h5 = yline(target_mod, ':', 'Model $0.5\,C(0)$', ...
        'Interpreter','latex', ...
        'LabelHorizontalAlignment','right');
else
    h5 = [];
end

% Empirical half-distance
if isfinite(Lhalf_emp)
    h6 = xline(Lhalf_emp, '--', ...
        sprintf('$\\xi_{1/2}^{(emp)} \\approx %.2f\\ \\mathrm{km}$', Lhalf_emp), ...
        'Interpreter','latex', ...
        'LabelOrientation','horizontal', ...
        'LabelVerticalAlignment','bottom', ...
        'LabelHorizontalAlignment','left');
else
    h6 = [];
end

% Model half-distance
if isfield(fit, 'Lhalf_mod') && isfinite(fit.Lhalf_mod)
    h7 = xline(fit.Lhalf_mod, ':', ...
        sprintf('$\\xi_{1/2}^{(mod)} \\approx %.2f\\ \\mathrm{km}$', fit.Lhalf_mod), ...
        'Interpreter','latex', ...
        'LabelOrientation','horizontal', ...
        'LabelVerticalAlignment','top', ...
        'LabelHorizontalAlignment','right');
else
    h7 = [];
end

xlabel('Lag distance $h$ (km)', 'Interpreter','latex');
ylabel('Covariance $C(h)$ (mGal$^2$)', 'Interpreter','latex');
title(['Empirical vs ', fit.model_name, ' covariance model'], 'Interpreter','latex');

% Legend handles / labels
legend_handles = [h1 h2 h3 h4];
legend_labels = { ...
    'Empirical covariance', ...
    sprintf('Valid bins ($N_{pairs} \\geq %d$)', min_pairs), ...
    [fit.model_name, ' fit'], ...
    '$0.5\,C(0)$'};

if ~isempty(h5)
    legend_handles = [legend_handles h5];
    legend_labels{end+1} = 'Model $0.5\,C(0)$';
end
if ~isempty(h6)
    legend_handles = [legend_handles h6];
    legend_labels{end+1} = '$\xi_{1/2}^{(emp)}$';
end
if ~isempty(h7)
    legend_handles = [legend_handles h7];
    legend_labels{end+1} = '$\xi_{1/2}^{(mod)}$';
end

legend(legend_handles, legend_labels, ...
    'Interpreter','latex', ...
    'Location','northeast');

% Annotation box
xl = xlim;
yl = ylim;

txt = { ...
    sprintf('Empirical $C(0)$: %.2f $\\mathrm{mGal}^2$', C0_emp), ...
    sprintf('Model $\\sigma^2$: %.2f $\\mathrm{mGal}^2$', fit.sigma2), ...
    sprintf('Model scale parameter: %.2f km', fit.a), ...
    sprintf('Weighted SSE: %.4g', fit.sse) ...
    };

if isfinite(Lhalf_emp)
    txt{end+1} = sprintf('Empirical $\\xi_{1/2}$: %.2f km', Lhalf_emp);
end
if isfield(fit, 'Lhalf_mod') && isfinite(fit.Lhalf_mod)
    txt{end+1} = sprintf('Model $\\xi_{1/2}$: %.2f km', fit.Lhalf_mod);
end

xl = xlim;
yl = ylim;

x_mid = mean(xl);
y_top = yl(1) + 0.85*(yl(2)-yl(1));

text(x_mid, y_top, txt, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'Interpreter','latex', ...
    'BackgroundColor','w');

hold off;
end