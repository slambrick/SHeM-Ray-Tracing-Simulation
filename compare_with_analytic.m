% Copyright (c) 2018, 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% Script for comparing the results from the ray-tracing method with the
% 'analytic' method of scattering off of randomly rough surfaces with an
% intrinsic specular scattering distribution. The agreement should be perfect at
% very low roughness, good at small roughness becoming increasingly poor as the
% roughness increases.

%% Parameters

% The range of roughness ratios to use.
ratios = 0.001:0.001:0.01;

% The number of surface elements to use
n_elements = 10000;

% The number of rays to use in the ray tracing model
n_rays = n_elements*2;

% The incidence angle
init_angle = -45;

% Should the mex file be recompiled
recompile = true;

% File to save results to
data_file = '../results1D/comparison_with_analytic_dfile001.mat';

%% Perform both calculations

loadpath

if recompile
   mex  CFLAGS='$CFLAGS -I mtwister -Wall' mexFiles/scatterRaysMex2D.c mexFiles/scattering2D.c ...
       mexFiles/scattering_processes2D.c mexFiles/small_functions2D.c ...
       mtwister/mtwister.c ...
       mexFiles/common_helpers.c
end

% Variables to store results in
analytic_thetas = zeros(length(ratios), n_elements);
analytic_weights = zeros(length(ratios), n_elements);
simulate_thetas = zeros(length(ratios), n_rays);

if isempty(gcp('nocreate'))
    parpool 
end
ppm = ParforProgressbar(length(ratios));
parfor i_=1:length(ratios)
    % Generate a surface
    [hs, xs] = rsgene1D(n_elements+1, (n_elements+1)/100, ratios(i_), 1);
    
    % Calculate the scattering distribution using the analtic method
    [analytic_thetas(i_,:), analytic_weights(i_,:)] = ...
        tracing2D.alalytic_random_scatter(xs, hs, init_angle, false);
    
    % Simulate the scattering distribution using the ray tracing method
    [simulate_thetas(i_,:), ~] = tracing2D.random_scatter('hs', hs, 'xs', xs, 'n_rays', ...
        n_rays, 'init_angle', -init_angle, 'scattering', 'specular');
    
    ppm.increment();
end
delete(ppm)

%% Plot results

% Loop through each ratio and plot both as two subplots on the same figure
for i_=1:length(ratios)
    figure
    subplot(2, 1, 1)
    weighted_histogram(analytic_thetas(i_,:), analytic_weights(i_,:), 30, 'pdf', false);
    xlabel('\theta')
    ylabel('P(\theta)')
    title('Analytic')
    xlim([-90 90])
    grid on

    subplot(2, 1, 2)
    histogram(simulate_thetas(i_,:), 30, 'normalization', 'pdf')
    xlabel('\theta')
    ylabel('P(\theta)')
    title('Ray tracing')
    xlim([-90 90])
    grid on
    
    sgtitle(['RMS height/correlation length = ' num2str(ratios(i_))])
end

% Save the results
save(data_file)

