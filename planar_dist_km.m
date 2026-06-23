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


function D = planar_dist_km(varargin)
%PLANAR_DIST_KM Compute local planar distance in km

deg2km = 111.195;

% Case 1: two matrices [lat lon]
if nargin == 2
    Z1 = varargin{1};
    Z2 = varargin{2};

    lat1 = Z1(:,1);
    lon1 = Z1(:,2);

    lat2 = Z2(:,1);
    lon2 = Z2(:,2);

    n1 = size(Z1,1);
    n2 = size(Z2,1);

    ref_lat = 0.5 * (mean(lat1) + mean(lat2));
    lon_scale = cosd(ref_lat) * deg2km;

    D = zeros(n1, n2);

    for i = 1:n1
        dx = (lon2 - lon1(i)) * lon_scale;
        dy = (lat2 - lat1(i)) * deg2km;

        D(i,:) = hypot(dx, dy);
    end

    % Case 2: explicit lat/lon inputs

elseif nargin == 4
    lat1 = varargin{1};
    lon1 = varargin{2};
    lat2 = varargin{3};
    lon2 = varargin{4};

    ref_lat = 0.5 * (mean(lat1(:)) + mean(lat2(:)));

    dx = (lon2 - lon1) .* cosd(ref_lat) * deg2km;
    dy = (lat2 - lat1) * deg2km;

    D = hypot(dx, dy);

else
    error(['Invalid number of inputs. Use either ', ...
        'planar_dist_km(Z1, Z2) or planar_dist_km(lat1, lon1, lat2, lon2).']);
end
end