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


% This main script was crated only for testing the functions, before using
% them in the final app

% !!!!!!! IT WAS CREATED FOR TESTING THE FUNCTIONS BEFORE THE FINAL GUI APP

% ########## Introntaction ##########
% Σε αυτό το script γίνεται η κατάλληλη επεξεργασία των δεδομένων τα οποία
% είτε δημιουργούνται από διαδικασία η οποία ελέγχεται από τον χρήστη ορίζαντας
% τις παραμέτρους που χρειάζονται για την κατασκευή του simulated
% (synthetic) field, είτε εισάγονται ως αρχείο πραγματικών δεδομένων. Τα
% δεδομένα είναι residuals gravity anomalies. Το πρόγραμμα προσδιορίζει το
% εμπειρικό μοντέλο των δεδομένων και προσαρμόζει ύστερα από επιλογή του
%χρήστη το κατάλληλο αναλυτικό μοντέλο. Υστερα έχει την δυνατότητα ο
%} χρήστης να πραγματοποιησει...

clc; clear; close all

set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


% ------ SETTINGS ------


enable_Sfield = true;      % true  -> synthetic workflow available
                           % false -> file-based workflow

do_model_comparison   = true;   % compare analytical covariance models
do_synthetic_lsc_test = true;   % only meaningful for synthetic field
do_user_prediction    = false;  % prediction on user-supplied file

modelID = 'gaussian';     % 'exponential', 'gaussian', 'hirvonen'

bin_width = [];
max_dist  = 150;
reg_rel   = 1e-10;

% ------ DATA PREPARATION ------


if enable_Sfield
    % Synthetic field
    n_points = 3000;
    lat_min  = 36.5;
    lat_max  = 38.5;
    lon_min  = 24.0;
    lon_max  = 26.0;

    % Field settings
    seed        = 15;
    grid_dim    = 200;
    sigma_cells = 15; 
    signal_std  = 25;

    data = Synthetic_field(n_points, lat_min, lat_max, lon_min, lon_max, seed, grid_dim, sigma_cells, signal_std);

else
    % Observation file
    obs_filename = 'Data_test.txt';
    data = read_gravity_data(obs_filename);
end

% Assign variables
id   = data.id;
lat  = data.lat;
lon  = data.lon;
grav = data.grav;
stdg = data.std_grav;


% ------ EMPIRICAL COVARIANCE + MODEL COMPARISON ------

if do_model_comparison

    results = Empirical_Model(lat, lon, grav, bin_width, max_dist);

    % Optional plot of empirical covariance only
    % plot_empirical_covariance(results);

    % --- Fit models on full dataset
    fit_exp = fit_exponential_model(results);
    fit_gau = fit_gaussian_model(results);
    fit_hir = fit_hirvonen_model(results);

    % Forsberg model
    results_planar = Empirical_Model_Planar(lat, lon, grav, [], 150);
    fit_for = dynamic_model_fitting_test(results_planar, 'forsbergplanar');

    fprintf('\nModel comparison on full dataset:\n');
    fprintf('Exponential SSE = %.3f\n', fit_exp.sse);
    fprintf('Gaussian    SSE = %.3f\n', fit_gau.sse);
    fprintf('Hirvonen    SSE = %.3f\n', fit_hir.sse);

    % Plot the model you want to inspect
    switch lower(modelID)
        case 'exponential'
            plot_empirical_vs_analytical(results, fit_exp);
        case 'gaussian'
            plot_empirical_vs_analytical(results, fit_gau);
        case 'hirvonen'
            plot_empirical_vs_analytical(results, fit_hir);
        otherwise
            error('Unknown modelID: %s', modelID);
    end

end


% ------ SYNTHETIC LSC VALIDATION ------

if enable_Sfield && do_synthetic_lsc_test

    train_fraction = 0.8;

    val = validate_synthetic_lsc(data, modelID, train_fraction, bin_width, max_dist, reg_rel);

    fprintf('\nValidation using %s model:\n', val.fit.model_name);
    fprintf('RMSE (vs observed test values) = %.3f mGal\n', val.rmse_obs);
    fprintf('MAE  (vs observed test values) = %.3f mGal\n', val.mae_obs);
    fprintf('Bias (vs observed test values) = %.3f mGal\n', val.bias_obs);

    if isfield(data, 'grav_true') && ~isempty(data.grav_true)
        fprintf('RMSE (vs true synthetic field) = %.3f mGal\n', val.rmse_true);
        fprintf('MAE  (vs true synthetic field) = %.3f mGal\n', val.mae_true);
        fprintf('Bias (vs true synthetic field) = %.3f mGal\n', val.bias_true);
    end

end


% ------ USER PREDICTION WORKFLOW ------

if ~enable_Sfield && do_user_prediction

    % 1. Build empirical covariance from observation dataset
    results_obs = Empirical_Model(data.lat, data.lon, data.grav, bin_width, max_dist);

    % 2. Dynamic fitting for selected model
    fit_obs = fit_covariance_model(results_obs, modelID);

    % 3. Plot empirical vs selected analytical model
    plot_empirical_vs_analytical(results_obs, fit_obs);

    % 4. Read prediction points file [id lat lon height]
    pred_filename = 'Prediction_points.txt';
    pred_data = read_prediction_points(pred_filename);

    % 5. Run LSC prediction
    pred_out = run_lsc_prediction_user(data, pred_data, fit_obs, reg_rel);

    % 6. Show first predicted values
    nshow = min(10, numel(pred_out.id));
    disp('First predicted values [id, predicted anomaly, prediction std]:')
    disp([pred_out.id(1:nshow), pred_out.g_hat(1:nshow), pred_out.std_hat(1:nshow)])

    % 7. Export results
    export_lsc_predictions(pred_out, 'lsc_predictions.csv');

    fprintf('\nLSC prediction file exported: lsc_predictions.csv\n');
end



% ------ EXPORT_LSC_PREDICTIONS ------
function export_lsc_predictions(pred_out, filename)

% OUTPUT COLUMNS:
%   [id, lat, lon, height, predicted_gravity, prediction_std]

    M = [pred_out.id, ...
         pred_out.lat, ...
         pred_out.lon, ...
         pred_out.h, ...
         pred_out.g_hat, ...
         pred_out.std_hat];

    writematrix(M, filename);
end