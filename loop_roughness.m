% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the Sub-beam Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% A script for scattering rays off randomly rough 1D surfaces, it loops through
% a number of different ratios between the rms height and the correlation length
% and saves the results to a data file.
clear;
close all;
clc

%% Parameters

% Values of the ratio to use
ratios = 0:0.02:0.4;

% The code cannot cope with a ratio of 0, so we make it small but finite
ratios(1) = 0.001;

% Number of surface elements to use
n_elements = 10000;

% Number of rays to use
n_rays = n_elements*2;

% Incidence angle
init_angle = 45;

% Recompile the mex file?
recompile = false;

% File to save the data to
file_name = 'increasing_roughness011.mat';

%% Perform simulations

if ~exist('../results1D/', 'dir')
    mkdir('../results1D/')
end
file_name = ['../results1D/' file_name];

addpath(genpath('functions'), genpath('classes'), genpath('ParforProgMon'))

mexCompile(recompile);

% Variables to store results
thetas = zeros(length(ratios), n_rays);
num_scatters = zeros(size(thetas));

tic

% Loop through all the different roughnesses, parrallized for loop
if isempty(gcp('nocreate'))
    parpool 
end
ppm = ParforProgressbar(length(ratios));
parfor i_=1:length(ratios)
    [thetas(i_,:), num_scatters(i_,:)] = tracing2D.random_scatter('ratio', ratios(i_), ...
        'n_rays', n_rays, 'Nelements', n_elements, 'init_angle', init_angle, ...
        'scattering', 'specular');
    ppm.increment();
end
delete(ppm);

current_pool = gcp('nocreate');
delete(current_pool);

t = toc;

disp(['Simulation took ' num2str(t/60) 'mins to complete.'])

%% Save results

save(file_name)

