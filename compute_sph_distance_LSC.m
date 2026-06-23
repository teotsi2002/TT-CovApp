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


function D = compute_sph_distance_LSC(lat1, lon1, lat2, lon2)
% COMPUTE_SPH_DISTANCE_SETS Great-circle distances between two point sets

R = 6371; % Earth radius in km

lat1 = lat1(:);
lon1 = lon1(:);
lat2 = lat2(:);
lon2 = lon2(:);

rlat1 = deg2rad(lat1);
rlon1 = deg2rad(lon1);
rlat2 = deg2rad(lat2);
rlon2 = deg2rad(lon2);

dlat = rlat2.' - rlat1;
dlon = rlon2.' - rlon1;

a = sin(dlat/2).^2 + cos(rlat1).*cos(rlat2.').*sin(dlon/2).^2;
a = min(max(a,0),1);   % numerical safety

D = 2 * R * atan2(sqrt(a), sqrt(1-a));
end