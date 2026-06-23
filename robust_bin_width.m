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


function [bin_width, avgNN, medianNN] = robust_bin_width(D, max_dist)
% ROBUST_BIN_WIDTH Estimate robust bin width from nearest-neighbor distances.
%
% The bin width must not be too small, otherwise the first empirical
% covariance bins become noisy and may exceed C(0).

Dtmp = D;
Dtmp(Dtmp == 0) = NaN;

min_d = min(Dtmp, [], 2, 'omitnan');

avgNN    = mean(min_d, 'omitnan');
medianNN = median(min_d, 'omitnan');

if ~isfinite(avgNN) || avgNN <= 0
    avgNN = medianNN;
end

if ~isfinite(medianNN) || medianNN <= 0
    medianNN = avgNN;
end

% Candidate 1: based on nearest-neighbor distance
bw_nn = 2.0 * medianNN;

% Candidate 2: avoid too many bins
targetBins = 100;
bw_bins = max_dist / targetBins;

% Final bin width
bin_width = max([bw_nn, bw_bins, eps]);
end