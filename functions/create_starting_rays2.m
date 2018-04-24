function [ray_pos, ray_dir, n_rays] = create_starting_rays2(ray_sep, ...
        pinhole_c, pinhole_r, multipl, ...
        plot_starting_positions, thePath, alpha)
% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% create_starting_rays2(...)
% 
%  Generates a distribution of ray positions and directions across the pinhole
%  given the location and size of the pinhole. The model used to generate
%  the directions of the rays assume a quadratic distribution of angles
%  alpha*theta^2.
%      The rays are aranged in a square grid across the circular pinhole
%  and the number of rays at each point is determined by a Gaussian distribution
%  centered on the centre of the pinhole. The number of rays at each point
%  is given by multipl. The positions are generated in a different coordinate
%  system and then nedd to be roated and shifted to match the main system used
%  in the simulation.
%      Rays are automatically transformmed into the coordinate system used in
%  the main simulation.
%
%  INPUTS:
%   ray_sep       - the seperation of the rays in the grid pattern (mm)
%   pinhole_c     - the centre of the pinhole (three element vector)
%   pinhole_r     - the radius of the pinhole (mm)
%   multipl       - the number of raysa at the centre point of the pinhole
%   plot_starting_postions  - make plots of the starting postions of the rays
%                             true/false
%   thePath       - The path to save the figures to
%   alpha         - Constant in the predicted quadratic angular distribution
%                   (optional), this is provided in mrad^-2
%
%  OUTPUTS:
%   ray_pos - Nx3 array, each row contains the position of a ray
%   ray_dir - Nx3 array, each row contains the direction of a ray, directions
%             are normalized
%   n_rays  - the number of rays that have been generated
    
    if nargin == 6
        alpha = 69.715282455;
    elseif nargin ~= 7
        error('Wrong number of arguments for create_starting_rays2()')
    end
    
    %% Generate the starting positions
    % ray_pos2 stores the locations of the rays in a coordinate system centered
    % on the center of the pinhole itself, at pi/2 to the simulation basis.
    ray_pos2 = [];
    
    n_rays = 1;
    for i_ = -pinhole_r:ray_sep:pinhole_r
        for j_ = -pinhole_r:ray_sep:pinhole_r
            if ((i_^2 + j_^2) < pinhole_r^2)
                % If inside the radius of the pinhole add to the inital
                % positions of the rays. The y position is fixed by the
                % definition of the y-axis.
                ray_pos2(n_rays,:) = [i_, 0, j_]; %#ok<AGROW>
                n_rays = n_rays + 1;
            end
        end
    end
    
    % Project the ray positions into a the coordinate basis being used for the
    % simulation
    ray_pos = ray_pos2;
    ray_pos(:,1) = (ray_pos2(:,1) - ray_pos2(:,2))/sqrt(2);
    ray_pos(:,2) = (ray_pos2(:,1) + ray_pos2(:,2))/sqrt(2);
    ray_pos(:,3) = ray_pos2(:,3);
    
    %% Plot the starting positions
    if plot_starting_positions
        figure
        plot(ray_pos(:,1), ray_pos(:,3), '.', 'color', [0.8 0.2 0])
        xlabel('x/mm')
        ylabel('z/mm')
        axis('equal')
        xlim([-0.002 0.002])
        ylim([-0.003 0.003])
        saveas(gcf, [thePath '/starting_positions1'], 'epsc')
        
        figure
        plot(ray_pos(:,1), ray_pos(:,2), '.', 'color', [0.8 0.2 0])
        hold on
        plot([-pinhole_r*sqrt(2), pinhole_r*sqrt(2)], [0, 0], 'color', [0.1 0.7 0.2], 'linewidth', 2)
        xlabel('x/mm')
        ylabel('y/mm')
        axis('equal')
        xlim([-0.004 0.004])
        ylim([-0.002 0.002])
        saveas(gcf, [thePath '/starting_positions2'], 'epsc')
        hold off
        
        figure
        plot(ray_pos(:,2), ray_pos(:,3), '.', 'color', [0.8 0.2 0])
        hold on
        plot([0,0], [-pinhole_r, pinhole_r], 'color', [0.1 0.7 0.2], 'linewidth', 3)
        xlabel('y/mm')
        ylabel('z/mm')
        axis('equal')
        xlim([-0.002 0.002])
        ylim([-0.003 0.003])
        saveas(gcf, [thePath '/starting_positions3'], 'epsc')
    end
    
    %% Move starting positions
    % Move the starting positions to the location of the pinhole
    ray_pos = bsxfun(@plus, ray_pos, pinhole_c);
    
    %% Repeat the starting positions
    ray_pos = repmat(ray_pos, multipl, 1);
    n_rays = size(ray_pos, 1);
    
    %% Generate the directions according to the predicted distribution of angles
    % random_dir uses y = alpha*x^2
    ray_dir2 = random_dir(n_rays, alpha);
    
    % Convert the directions into the correct coordinate system
    ray_dir = ray_dir2;
    ray_dir(:,1) = (ray_dir2(:,1) - ray_dir2(:,2))/sqrt(2);
    ray_dir(:,2) = (ray_dir2(:,1) + ray_dir2(:,2))/sqrt(2);
    ray_dir(:,3) = ray_dir2(:,3);
end







