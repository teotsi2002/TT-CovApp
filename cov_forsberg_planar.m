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


function C = cov_forsberg_planar(d_km, h1_km, h2_km, C0, D, T)
% COV_FORSBERG_PLANAR Forsberg planar covariance model for gravity anomalies

% ------ Basic parameter checks ------
if C0 <= 0 || D <= 0 || T <= 0
    error('Forsberg parameters must satisfy C0 > 0, D > 0, T > 0.');
end

% ------ Auxiliary depth terms ------
depth_terms = [D, D + T, D + 2*T, D + 3*T];
coeffs      = [1, -3, 3, -1];

if any(depth_terms <= 0) || any(~isfinite(depth_terms))
    error('Invalid Forsberg parameters: depth terms must be positive and finite.');
end

% ------ Normalization factor so that C(0,0,0) = C0 ------
norm_factor = C0 / log((depth_terms(2)^3 * depth_terms(4)) / ...
    (depth_terms(1)   * depth_terms(3)^3));

% ------ Heights enter through their sum ------
h_sum = h1_km + h2_km;

% ------ Covariance summation ------
C_sum = zeros(size(d_km));

for k = 1:4
    z_k = h_sum + depth_terms(k);
    r_k = sqrt(d_km.^2 + z_k.^2);

    C_sum = C_sum + coeffs(k) * (-log(z_k + r_k));
end

% ------ Final covariance ------
C = norm_factor * C_sum;
end