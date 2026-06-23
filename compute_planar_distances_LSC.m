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


function D = compute_planar_distance_LSC(lat1, lon1, lat2, lon2)
%COMPUTE_PLANAR_DISTANCE_LSC Planar horizontal distances between two point sets

R = 6371.0; % km

lat1 = lat1(:);
lon1 = lon1(:);
lat2 = lat2(:);
lon2 = lon2(:);

lat0 = mean([lat1; lat2], 'omitnan');
km_per_deg_lat = 111.32;
km_per_deg_lon = 111.32 * cosd(lat0);

x1 = lon1 * km_per_deg_lon;
y1 = lat1 * km_per_deg_lat;

x2 = lon2 * km_per_deg_lon;
y2 = lat2 * km_per_deg_lat;

dx = x2.' - x1;
dy = y2.' - y1;

D = hypot(dx, dy);
end