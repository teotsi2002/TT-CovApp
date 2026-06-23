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


function results = Empirical_Model(lat, lon, grav, bin_width, max_dist)
% EMPIRICAL_MODEL Empirical covariance from spherical distances.

% ------ Input formatting ------
lat  = lat(:);
lon  = lon(:);
grav = grav(:);

N = numel(grav);

if numel(lat) ~= N || numel(lon) ~= N
    error('lat, lon, and grav must have the same length.');
end

if N < 3
    error('At least 3 points are required for empirical covariance estimation.');
end

% ------ Remove invalid observations ------
validObs = isfinite(lat) & isfinite(lon) & isfinite(grav);

lat  = lat(validObs);
lon  = lon(validObs);
grav = grav(validObs);

N = numel(grav);

if N < 3
    error('Not enough valid observations after removing NaN/Inf values.');
end

% ------ Center gravity values ------
grav = grav - mean(grav, 'omitnan');

% ------ Zero-lag variance ------
% C(0) is estimated separately and is not taken from the first distance bin.
C0_emp = mean(grav.^2, 'omitnan');

% ------ Distance matrix ------
D = compute_sph_distance(lat, lon);

% ------ Unique pairs only ------
idx = triu(true(N), 1);

d_pairs = D(idx);

prod_matrix = grav * grav.';
prod_pairs = prod_matrix(idx);

valid_pairs = isfinite(d_pairs) & isfinite(prod_pairs);

d_pairs    = d_pairs(valid_pairs);
prod_pairs = prod_pairs(valid_pairs);

if isempty(d_pairs)
    error('No valid point pairs available for empirical covariance estimation.');
end

% ------ Maximum lag distance ------
if nargin < 5 || isempty(max_dist)
    max_dist = prctile(d_pairs, 99);
end

if ~isfinite(max_dist) || max_dist <= 0
    error('max_dist must be positive and finite.');
end

keep = d_pairs <= max_dist;

d_pairs    = d_pairs(keep);
prod_pairs = prod_pairs(keep);

if isempty(d_pairs)
    error('No point pairs remain after applying max_dist.');
end

% ------ Automatic bin width ------
if nargin < 4 || isempty(bin_width)
    [bin_width, avgNN, medianNN] = robust_bin_width(D, max_dist);
else
    avgNN = NaN;
    medianNN = NaN;
end

if ~isfinite(bin_width) || bin_width <= 0
    error('bin_width must be positive and finite.');
end

% Avoid an excessive number of too narrow bins
maxBinsAllowed = 120;

if ceil(max_dist / bin_width) > maxBinsAllowed
    bin_width = max_dist / maxBinsAllowed;
end

% ------ Define bins ------
edges = 0:bin_width:(max_dist + bin_width);

nBins = numel(edges) - 1;

bin_centers = zeros(nBins,1);
mean_dist   = nan(nBins,1);
emp_cov     = nan(nBins,1);
pair_count  = zeros(nBins,1);

for k = 1:nBins

    if k < nBins
        in_bin = d_pairs >= edges(k) & d_pairs < edges(k+1);
    else
        in_bin = d_pairs >= edges(k) & d_pairs <= edges(k+1);
    end

    pair_count(k) = sum(in_bin);
    bin_centers(k) = 0.5 * (edges(k) + edges(k+1));

    if pair_count(k) > 0
        mean_dist(k) = mean(d_pairs(in_bin), 'omitnan');
        emp_cov(k)   = mean(prod_pairs(in_bin), 'omitnan');
    end
end

% ------ Minimum pairs per bin ------
% 30 is too small for large real datasets. Use a dynamic threshold,
% but keep it moderate so that the first part of the curve is not lost.
totalPairsUsed = numel(d_pairs);

min_pairs_bin = max(30, min(300, round(0.0005 * totalPairsUsed)));

valid_bins = isfinite(mean_dist) & isfinite(emp_cov) & ...
             pair_count >= min_pairs_bin;

% ------ Main empirical sequence used for plotting/fitting ------
lags = [0; mean_dist(:)];
C_emp = [C0_emp; emp_cov(:)];
pair_count_with_zero = [N; pair_count(:)];

lags_valid_with_zero = [0; mean_dist(valid_bins)];
C_emp_valid_with_zero = [C0_emp; emp_cov(valid_bins)];
pair_count_valid_with_zero = [N; pair_count(valid_bins)];

% ------ Empirical half-covariance distance ------
target = 0.5 * C0_emp;
Lhalf_emp = NaN;

x = lags_valid_with_zero;
y = C_emp_valid_with_zero;

finiteXY = isfinite(x) & isfinite(y);

x = x(finiteXY);
y = y(finiteXY);

if numel(x) >= 2

    % Start from h = 0 and find first crossing of 0.5*C(0)
    k = find(y <= target, 1, 'first');

    if ~isempty(k) && k > 1
        x1 = x(k-1); x2 = x(k);
        y1 = y(k-1); y2 = y(k);

        if y2 ~= y1
            Lhalf_emp = x1 + (target - y1) * (x2 - x1) / (y2 - y1);
        else
            Lhalf_emp = x2;
        end
    end
end

% ------ Output ------
results.D = D;
results.d_pairs = d_pairs;
results.prod_pairs = prod_pairs;

results.bin_width = bin_width;
results.avgNN = avgNN;
results.medianNN = medianNN;
results.max_dist = max_dist;

results.edges = edges;
results.bin_centers = bin_centers;
results.mean_dist = mean_dist;
results.emp_cov = emp_cov;
results.pair_count = pair_count;

results.valid_bins = valid_bins;

% Valid bins without zero-lag, for compatibility with existing fit functions
results.lags_valid = mean_dist(valid_bins);
results.C_emp_valid = emp_cov(valid_bins);
results.pair_count_valid = pair_count(valid_bins);

% Valid bins with zero-lag, for plotting and diagnostics
results.lags_valid_with_zero = lags_valid_with_zero;
results.C_emp_valid_with_zero = C_emp_valid_with_zero;
results.pair_count_valid_with_zero = pair_count_valid_with_zero;

results.min_pairs_bin = min_pairs_bin;
results.nBins = nBins;

results.C0_emp = C0_emp;
results.lags = lags;
results.C_emp = C_emp;
results.pair_count_with_zero = pair_count_with_zero;

results.Lhalf_emp = Lhalf_emp;
end