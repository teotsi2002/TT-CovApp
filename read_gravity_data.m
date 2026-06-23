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


function data = read_gravity_data(filename)
%READ_GRAVITY_DATA Read gravity anomaly data from file.

% Accepted formats:
% [id, latitude, longitude, gravity_anomaly]
% [id, latitude, longitude, gravity_anomaly, std_gravity]

% The function first tries readmatrix() for purely numeric files.
% If that is not suitable, it falls back to readtable() for files with headers.

raw = [];

% ------ First attempt: numeric matrix ------
try
    raw = readmatrix(filename);

    % Remove completely empty rows
    raw = raw(~all(isnan(raw),2), :);

    % Check if readmatrix result is usable
    if isempty(raw)
        raw = [];
    else
        [~, nCols_tmp] = size(raw);

        if ~(nCols_tmp == 4 || nCols_tmp == 5)
            raw = [];
        elseif any(any(isnan(raw(:,1:min(4,nCols_tmp)))))
            % Probably a file with headers or mixed content
            raw = [];
        end
    end

catch
    raw = [];
end

% ------ Second attempt: table (for files with headers) ------
if isempty(raw)
    try
        T = readtable(filename);
        raw = table2array(T);

        % Remove completely empty rows
        raw = raw(~all(isnan(raw),2), :);

        if isempty(raw)
            error('The file is empty.');
        end

    catch ME
        error('Unable to read file "%s". MATLAB message: %s', filename, ME.message);
    end
end

% ------ Format check ------
[nRows, nCols] = size(raw);

if ~(nCols == 4 || nCols == 5)
    error(['Invalid file format. Accepted formats are only: ', ...
        '[id, latitude, longitude, gravity_anomaly] or ', ...
        '[id, latitude, longitude, gravity_anomaly, std_gravity].']);
end

% First 4 columns must exist and be numeric
required_part = raw(:,1:4);

if any(any(isnan(required_part)))
    error(['The first 4 required columns contain missing or non-numeric values. ', ...
        'Please check the file contents.']);
end

% ------ Assign variables ------
data.id   = raw(:,1);
data.lat  = raw(:,2);
data.lon  = raw(:,3);
data.grav = raw(:,4);

if nCols == 5
    data.std_grav = raw(:,5);
else
    data.std_grav = [];
end

data.n_points = nRows;
data.n_cols   = nCols;

% ------ Value checks ------

% ID must be integer
if any(mod(data.id,1) ~= 0)
    error('ID values must be integers.');
end

% No duplicate IDs
if numel(unique(data.id)) < numel(data.id)
    error('Duplicate IDs were found in the file.');
end

% Latitude range
if any(data.lat < -90 | data.lat > 90)
    error('Latitude values must be in the range [-90, 90] degrees.');
end

% Longitude range
if any(data.lon < -180 | data.lon > 180)
    error('Longitude values must be in the range [-180, 180] degrees.');
end

% Gravity anomaly values
if any(isnan(data.grav))
    error('Gravity anomaly values contain invalid entries.');
end

% Standard deviation values
if nCols == 5
    if any(isnan(data.std_grav))
        error('Standard deviation values contain invalid entries.');
    end

    if any(data.std_grav <= 0)
        error('Standard deviation values must be positive.');
    end
end

% No duplicate coordinates
coords = [data.lat, data.lon];
if size(unique(coords, 'rows'), 1) < size(coords, 1)
    error('Duplicate coordinate pairs (latitude, longitude) were found.');
end

end