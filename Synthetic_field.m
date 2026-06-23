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


function data = Synthetic_field(n_points, lat_min, lat_max, lon_min, lon_max, seed, grid_dim, sigma_cells, signal_std)
%SYNTHETIC_FIELD Generate synthetic gravity anomaly data with spatial correlation

rng(seed);

% ------ Region boundaries ------
buffer   = 0.5;

% Buffered grid
lon_grid   = linspace(lon_min - buffer, lon_max + buffer, grid_dim);
lat_grid   = linspace(lat_min - buffer, lat_max + buffer, grid_dim);
[Lon, Lat] = meshgrid(lon_grid, lat_grid);


% ------ Random correlated signal ------
random_noise = randn(size(Lon));

filter_half_size = ceil(3 * sigma_cells);

x_filter = -filter_half_size : filter_half_size;
y_filter = -filter_half_size : filter_half_size;
[X, Y]   = meshgrid(x_filter, y_filter);

gauss_filter = exp(-(X.^2 + Y.^2) / (2 * sigma_cells^2));
gauss_filter = gauss_filter / sum(gauss_filter(:));

% Convolution
signal_grid = conv2(random_noise, gauss_filter, 'same');

% Scale to mGal
signal_grid = (signal_grid - mean(signal_grid(:))) / std(signal_grid(:));
signal_grid = signal_std * signal_grid;

% ------ Sample random observation points ------

lat_points = lat_min + (lat_max - lat_min) * rand(n_points,1);
lon_points = lon_min + (lon_max - lon_min) * rand(n_points,1);

grav_true = interp2(Lon, Lat, signal_grid, lon_points, lat_points, 'cubic');

% ------ Add measurement noise ------

sigma_grav  = 1 + 1.5 * rand(n_points,1);
white_noise = randn(n_points,1) .* sigma_grav;

grav_obs = grav_true + white_noise;
grav_obs = grav_obs - mean(grav_obs);

%  ------ Output ------
data.id        = (1:n_points)';
data.lat       = lat_points;
data.lon       = lon_points;
data.grav      = grav_obs;
data.std_grav  = sigma_grav;
end