% Copyright (c) 2020, Dan Seremet.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.
% composite.m
%
% This script reads the data outputs of three simulations of the same sample,
% a 1x1 square composed of two halves:
%   1.  using the material-sensitive code, half of the sample was assigned
%       a cosine distribution, and the other a uniform distribution.
%   2.  using the old original code, the sample was simulated
%       with a uniform distribution
%   3.  the same for a cosine distirbution
%
% The goal is to create a composite of the outputs of simulations 2 and 3 and
% show how similar it is to that of simulation 1.
%

sims_root = '../results';

% load in the simulation data
sim_material_code = load(...
    fullfile(sims_root, 'two_halves/scatteringData.mat')).simulationData;
sim_uniform = load(...
    fullfile(sims_root, 'two_halves_old_uniform/scatteringData.mat')).simulationData;
sim_cosine = load(...
    fullfile(sims_root, 'two_halves_old_cosine/scatteringData.mat')).simulationData;

% take the average counters of the left and right halves of the material-sensitive code
left_half = slice_by_coords(sim_material_code.cntrSum{1},...
    [-0.55, 0.55], [-0.55, 0.55], [-.5, -0.01], [-.5, .5]);

right_half = slice_by_coords(sim_material_code.cntrSum{1},...
    [-0.55, 0.55], [-0.55, 0.55], [0.01, 0.5], [-.5, .5]);

disp("=== average brightness for new code ===")
avg_left = mean(left_half, 'all');
std_left = std(left_half, 0, 'all');
avg_right = mean(right_half, 'all');
std_right = std(right_half, 0, 'all');
disp(['uniform (left): ' num2str(avg_left) ' +/- ' num2str(std_left)]);
disp(['cosine (right): ' num2str(avg_right) ' +/- ' num2str(std_right)]);

% take the left half of the uniform and right half of the cosine sample
uniform_cntr_sum = sim_uniform.cntrSum{1}(:, 1:111);
cosine_cntr_sum = sim_cosine.cntrSum{1}(:, 112:221);
cntr_sums = { cat(2, uniform_cntr_sum, cosine_cntr_sum) };

% slice out the empty space so we average brightness only over sample pixels
uniform_cntr_sample = slice_by_coords(cntr_sums{1},...
    [-0.55, 0.55], [-0.55, 0.55], [-.5, -0.01], [-.5, .5]);

cosine_cntr_sample = slice_by_coords(cntr_sums{1},...
    [-0.55, 0.55], [-0.55, 0.55], [0.01, 0.5], [-.5, .5]);

disp("=== average brightness for composite of old code ===")
avg_left = mean(uniform_cntr_sample, 'all');
std_left = std(uniform_cntr_sample, 0, 'all');
avg_right = mean(cosine_cntr_sample, 'all');
std_right = std(cosine_cntr_sample, 0, 'all');
disp(['uniform (left): ' num2str(avg_left) ' +/- ' num2str(std_left)]);
disp(['cosine (right): ' num2str(avg_right) ' +/- ' num2str(std_right)]);

% the max number of scatters is the number of counters at each pixel
max_scatters = size(sim_uniform.counters{1}, 1);

% construct the new composite simulationData object to display a figure
counters = { cat(3, sim_uniform.counters{1}(:,:,1:111), sim_cosine.counters{1}(:,:,112:221)) };

composite_data = sim_uniform;
composite_data.counters = counters;
composite_data.cntrSum = cntr_sums;

composite_data.imageAll();
sim_material_code.imageAll();

% composite_data = RectangleInfo(...
%     counters, sim_uniform.num_killed, sim_uniform.sample_surface, sim_uniform.xrange,...
%     sim_uniform.zrange, sim_uniform.raster_movment_x, sim_uniform.raster_movment_z,...
%     sim_uniform.rays_per_pixel, sim_uniform.n_effuse, sim_uniform.time,...
%     sim_uniform.time_estimate, sim_uniform.counter_effusive, sim_uniform.n_detector,...
%     max_scatters, sim_uniform.dist_to_sample, sim_uniform.beam_param);
