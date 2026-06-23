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


function results = Empirical_Model_Planar(lat, lon, grav, bin_width, max_dist)
% EMPIRICAL_MODEL_PLANAR Empirical covariance using planar Euclidean distances

% Input formatting
lat  = lat(:);
lon  = lon(:);
grav = grav(:);

N = numel(grav);

if numel(lat) ~= N || numel(lon) ~= N
    error('lat, lon, and grav must have the same length.');
end

% Center gravity values
grav = grav - mean(grav, 'omitnan');

% Planar distance matrix
D = planar_dist_km(lat, lon, lat.', lon.');

% Unique pairs only
idx = triu(true(N), 1);

d_pairs = D(idx);

prod_matrix = grav * grav.';
prod_pairs  = prod_matrix(idx);

valid_pairs = ~isnan(d_pairs) & ~isnan(prod_pairs);
d_pairs = d_pairs(valid_pairs);
prod_pairs = prod_pairs(valid_pairs);

% Automatic bin width
if nargin < 4 || isempty(bin_width)
    [bin_width, avgNN] = robust_bin_width(D);
else
    avgNN = NaN;
end

% Maximum lag distance
if nargin < 5 || isempty(max_dist)
    max_dist = prctile(d_pairs, 99);
end

keep = d_pairs <= max_dist;
d_pairs = d_pairs(keep);
prod_pairs = prod_pairs(keep);

% Define bins
edges = 0:bin_width:(max_dist + bin_width);
nBins = numel(edges) - 1;

bin_centers = zeros(nBins,1);
mean_dist   = nan(nBins,1);
emp_cov     = nan(nBins,1);
pair_count  = zeros(nBins,1);

for k = 1:nBins
    in_bin = d_pairs >= edges(k) & d_pairs < edges(k+1);

    pair_count(k) = sum(in_bin);
    bin_centers(k) = 0.5 * (edges(k) + edges(k+1));

    if pair_count(k) > 0
        mean_dist(k) = mean(d_pairs(in_bin));
        emp_cov(k)   = mean(prod_pairs(in_bin));
    end
end

% Zero-lag variance
C0_emp = var(grav, 1);

lags = [0; mean_dist(:)];
C_emp = [C0_emp; emp_cov(:)];
pair_count_with_zero = [N; pair_count(:)];

% Valid bins
min_pairs_bin = 30;
valid_bins = ~isnan(emp_cov) & (pair_count >= min_pairs_bin);

% Empirical half-covariance distance
target = 0.5 * C0_emp;
Lhalf_emp = NaN;

if any(valid_bins)
    x = mean_dist(valid_bins);
    y = emp_cov(valid_bins);

    k = find(y <= target, 1, 'first');

    if ~isempty(k) && k > 1
        x1 = x(k-1); x2 = x(k);
        y1 = y(k-1); y2 = y(k);

        if y2 ~= y1
            Lhalf_emp = x1 + (target - y1) * (x2 - x1) / (y2 - y1);
        end
    end
end

% ------ Output ------
results.D = D;
results.d_pairs = d_pairs;
results.prod_pairs = prod_pairs;

results.bin_width = bin_width;
results.avgNN = avgNN;
results.max_dist = max_dist;
results.nBins = nBins;

results.edges = edges;
results.bin_centers = bin_centers;
results.mean_dist = mean_dist;
results.emp_cov = emp_cov;
results.pair_count = pair_count;

results.valid_bins = valid_bins;
results.lags_valid = mean_dist(valid_bins);
results.C_emp_valid = emp_cov(valid_bins);
results.pair_count_valid = pair_count(valid_bins);
results.min_pairs_bin = min_pairs_bin;

results.C0_emp = C0_emp;
results.lags = lags;
results.C_emp = C_emp;
results.pair_count_with_zero = pair_count_with_zero;

results.Lhalf_emp = Lhalf_emp;
results.distance_type = 'planar';
end