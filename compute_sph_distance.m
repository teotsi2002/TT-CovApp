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


function sph_D = compute_sph_distance(lat, lon)
% COMPUTE_SPHERICAL_DISTANCES Compute great-circle distances (km)

R = 6371; % Earth radius in km

% Ensure column vectors
lat = lat(:);
lon = lon(:);

% Convert to radians
rlat1 = deg2rad(lat);
rlon1 = deg2rad(lon);

% Create pairwise matrices using transpose
rlat2 = rlat1.';
rlon2 = rlon1.';

dlat = rlat2 - rlat1;
dlon = rlon2 - rlon1;

a = sin(dlat/2).^2 + cos(rlat1).*cos(rlat2).*sin(dlon/2).^2;
a = min(max(a,0),1); % numerical safety

sph_D = 2 * R * atan2(sqrt(a), sqrt(1-a));
end