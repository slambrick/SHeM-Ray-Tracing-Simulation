% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
%
% This file is for testing how the simulated pinhole beam spreads. It
% constructs, using the function create_starting_rays2() and random_dir(), rays
% that follow a predicted angular distribution that assumes the skimmer to be
% the source.
%
% The distribution is proportional to alpha*theta^2.
clear;

%% Path to files
addpath('functions', 'classes');

%% Parameters
ray_sep = 0.000002;    % Seperation of rays in the pinhole
pinhole_c = [0 0 0];   % Centre of the pinhole
pinhole_r = 0.00025;   % Radius of the pinhole
multipl = 1;           % Number of rays at each point in the pinhole
distance = 2.121;          % The distance (mm) to propogate the beam
alpha = 69.715282455;  % The constant in the angular distribution alpha*theta^2
save_to_file = false;  % Should the data be saved to a text file?
make_figure = true;    % Should a figure be made?
save_figure = false;   % should the figure be saved?

% If this is true the rays are parrallel and are not given directions according
% to the predicted distribution, hence they just follow straight lines.
parallel = false;

% Directory to save the figures and data to.
thePath = 'figures/beam_spread/beamSpread2';
texFname = 'beam_spread_positions2.csv';


%% Generate data
% Create starting rays, create_starting_rays produces rays that all have the
% same direction while create_starting_rays2 gives them a direction according to
% the predicted distribution
if parallel
    [ray_pos, ray_dir, n_rays] = create_starting_rays(ray_sep, pinhole_c, ...
        pinhole_r, multipl, false, thePath); %#ok<*UNRCH>
else    
    [ray_pos, ray_dir, n_rays] = create_starting_rays2(ray_sep, pinhole_c, ...
        pinhole_r, multipl, false, thePath, alpha);
end

% Need to transform the direction and position into another basis to make it 
% easier to plot
ray_pos2(:,1) = (ray_pos(:,1) + ray_pos(:,2))/sqrt(2);
ray_pos2(:,2) = (ray_pos(:,1) - ray_pos(:,2))/sqrt(2);
ray_pos2(:,3) = ray_pos(:,3);
ray_pos = ray_pos2;

ray_dir2(:,1) = (ray_dir(:,1) + ray_dir(:,2))/sqrt(2);
ray_dir2(:,2) = (ray_dir(:,1) - ray_dir(:,2))/sqrt(2);
ray_dir2(:,3) = ray_dir(:,3);
ray_dir = ray_dir2;

% Propogate the rays
new_pos = ray_pos + distance*ray_dir;

%% Plot/output data
% Plot the position of the rays
if ~exist(thePath, 'dir')
    mkdir(thePath)
end
if make_figure
    plot_pos = new_pos;
    figure
    plot(plot_pos(:,1), plot_pos(:,3), '.', 'Color', 'red', 'MarkerSize', 4);
    hold on 
    t = linspace(0,2*pi,100);
    plot(pinhole_r*sin(t),pinhole_r*cos(t), 'Color', 'blue', 'LineWidth', 3);
    xlabel('x/mm');
    ylabel('y/mm');
    axis('equal');
    set(gcf, 'PaperPosition', [0 0 12 11]);
    set(gcf, 'PaperSize', [12 11])
    if save_figure
        saveas(gcf, [thePath '/ray_positions'], 'epsc')
    end
end

% Output data of the final positions
if save_to_file
    dlmwrite([thePath '/' texFname], new_pos)
end

