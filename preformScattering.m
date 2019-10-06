% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

clear

%% Parameters

% Starting direction of the rays
init_theta = 45;
init_dir = [sind(init_theta), -cosd(init_theta), 0];

% Number of rays to trace
n_rays = 10000000;

% The maximum number of scattering events a ray is allowed to undergo
maxScatter = 1000;

% Specify the sample to scatter off and manipulate so you get the right
% spot.
sample_fname = 'deep_trench_sample.stl';

% Should the mex files be recompiled
recompile = false;

%% Perform calculation

% Paths to functions
addpath('stlread', 'functions', 'functions/interface_functions', 'classes', ...
        'mexFiles', 'DylanMuir-ParforProgMon-a80c9e9', 'functions/standard_samples', ...
        'surf2stl');

% Manipulate the sample and beam starting point
sample_surface = inputSample('fname', sample_fname, 'dontMeddle', true, ...
    'scattering', 1);

% Choose a spot on the sample then displace backwards
init_pos = [0, 0.5, 0];
init_pos = init_pos - 10*init_dir;

if recompile
    mexCompile();
end

% Main computation
[killed, numScattersRay, final_pos, final_dir] = ...
    distributionCalc('sample_surface', sample_surface, 'maxScatter', maxScatter, ...
                     'nrays', n_rays, 'start_pos', init_pos, 'start_dir', init_dir);

%% Analyse results

% Get the outgoing directions
ind = numScattersRay > 1;
[az, el, r] = cart2sph(final_dir(1,ind), final_dir(3,ind), final_dir(2,ind));

% Select the multiply scattered outgoing directions
th = (pi/2 - el)*180/pi;
[n, c] = hist3([th', az'*180/pi + 180], {2.5:5:87.5, 2.5:5:357.5}, 'Normalization', 'pdf');
n2 = zeros(size(n));
for i_=1:length(c{1})
    if c{1}(i_) ~= 0
        n2(i_,:) = n(i_,:)/sind(c{1}(i_));
    else
        n2(i_,:) = 0;
    end
end
n2 = n2/max(n2(:));
figure
polarPcolor(90*(c{1} - 2.5)/max(c{1} - 2.5), 360*(c{2} - 2.5)/max(c{2} - 2.5), n2', 'colBar', false)

% Polar density plot from the function
th = 2.5:5:87.5;
phi = -177.5:5:177.5;
th = th*pi/180;
phi = phi*pi/180;
[ths, phis] = meshgrid(th, phi);
dens = cosd(ths);
polarPcolor((th*180/pi - 2.5)*90/max(th*180/pi - 2.5), (phi*180/pi + 180 - 2.5)*360/max(phi*180/pi + 180 - 2.5), dens, 'colBar', false)



