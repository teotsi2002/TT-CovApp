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


function pred_data = read_prediction_points(filename)
%READ_PREDICTION_POINTS Read prediction point file

raw = [];

% ------ First attempt: readmatrix ------
try
    raw = readmatrix(filename);

    % Remove completely empty rows
    raw = raw(~all(isnan(raw),2), :);

    if isempty(raw)
        raw = [];
    end
catch
    raw = [];
end

% ------ Second attempt: readtable ------
if isempty(raw)
    try
        T = readtable(filename);
        raw = table2array(T);

        raw = raw(~all(isnan(raw),2), :);

        if isempty(raw)
            error('The prediction file is empty.');
        end
    catch ME
        error('Unable to read prediction file "%s". MATLAB message: %s', filename, ME.message);
    end
end

% ------ Format check ------
[nRows, nCols] = size(raw);

if nCols < 4
    error('Prediction file must contain at least 4 columns: [id, latitude, longitude, height].');
end

req = raw(:,1:4);

if any(any(isnan(req)))
    error('Prediction file contains missing or non-numeric values in the first 4 columns.');
end

% Assign fields
pred_data.id  = raw(:,1);
pred_data.lat = raw(:,2);
pred_data.lon = raw(:,3);
pred_data.h   = raw(:,4);

pred_data.n_points = nRows;

% ------Checks ------
if any(mod(pred_data.id,1) ~= 0)
    error('Prediction point IDs must be integers.');
end

if numel(unique(pred_data.id)) < numel(pred_data.id)
    error('Duplicate IDs were found in prediction file.');
end

if any(pred_data.lat < -90 | pred_data.lat > 90)
    error('Latitude values must be in the range [-90, 90] degrees.');
end

if any(pred_data.lon < -180 | pred_data.lon > 180)
    error('Longitude values must be in the range [-180, 180] degrees.');
end
end