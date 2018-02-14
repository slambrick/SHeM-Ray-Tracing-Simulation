% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

function [ray_pos, ray_dir, n_rays] = create_starting_rays(ray_sep, ...
        pinhole_c, pinhole_r, multipl, plot_starting_positions, ...
        thePath, FWHM_density, plot_density)
% create_starting_rays(...)
%
% Generates a distribution of ray positions and directions across the pinhole
% given the location and size of the pinhole. The rays are uniformally created
% by this function and are all created with the same initial direction.
%     The rays are aranged in a square grid across the circular pinhole
% and the number of rays at each point is determined by a Gaussian distribution
% centered on the centre of the pinhole. The number of rays at each point
% is given by multipl. The positions are generated in a different coordinate
% system and then nedd to be roated and shifted to match the main system used
% in the simulation.
% 
% INPUTS:
%  ray_sep       - the seperation of the rays in the grid pattern (mm)
%  pinhole_c     - the centre of the pinhole (three element vector)
%  pinhole_r     - the radius of the pinhole (mm)
%  multipl       - the number of raysa at the centre point of the pinhole
%  plot_starting_postions  - make a plot of the ray starting positions in the
%                            pinhole. true/false
%  thePath       - The path to save the figures to
%  FWHM_density  - The FWHM of the Gaussian beam density
%  plot_density  - Should a plot of the beam density be made, true/false
%
% OUTPUTS:
%  ray_pos - nx3 array, each row contains the position of a ray
%  ray_dir - nx3 array, each row contains the direction of a ray, directions
%            are normalized
%  n_rays  - the number of rays that have been generated

    
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
    ray_pos(:,1) = (ray_pos2(:,1) - ray_pos2(:,2))/sqrt(2);
    ray_pos(:,2) = (ray_pos2(:,1) + ray_pos2(:,2))/sqrt(2);
    ray_pos(:,3) = ray_pos2(:,3);
     
    if plot_starting_positions
        figure
        plot(ray_pos(:,1), ray_pos(:,3), '.', 'color', [0.8 0.2 0])
        xlabel('x/mm')
        ylabel('z/mm')
        axis('equal')
        xlim([-1.5*pinhole_r 1.5*pinhole_r])
        ylim([-1.5*pinhole_r 1.5*pinhole_r])
        saveas(gcf, [thePath '/starting_positions1'], 'epsc')
        
        figure
        plot(ray_pos(:,1), ray_pos(:,2), '.', 'color', [0.8 0.2 0])
        hold on
        plot([-pinhole_r*sqrt(2), pinhole_r*sqrt(2)], [0, 0], 'color', [0.1 0.7 0.2], 'linewidth', 2)
        xlabel('x/mm')
        ylabel('y/mm')
        axis('equal')
        xlim([-1.5*pinhole_r 1.5*pinhole_r])
        ylim([-1.5*pinhole_r 1.5*pinhole_r])
        saveas(gcf, [thePath '/starting_positions2'], 'epsc')
        hold off
        
        figure
        plot(ray_pos(:,2), ray_pos(:,3), '.', 'color', [0.8 0.2 0])
        hold on
        plot([0,0], [-pinhole_r, pinhole_r], 'color', [0.1 0.7 0.2], 'linewidth', 3)
        xlabel('y/mm')
        ylabel('z/mm')
        axis('equal')
        xlim([-1.5*pinhole_r 1.5*pinhole_r])
        ylim([-1.5*pinhole_r 1.5*pinhole_r])
        saveas(gcf, [thePath '/starting_positions3'], 'epsc')
    end
    
    ray_pos = bsxfun(@plus, ray_pos, pinhole_c);
    
    % Distribution of position by a Gaussian distribution
    sigma = FWHM_density/(2*sqrt(2*log(2)));
    intensity = mvnpdf([ray_pos2(:,1), ray_pos2(:,3)], [0 0], [sigma^2 0; 0 sigma^2]);
    intensity = intensity/max(intensity);
    nums = round(intensity*multipl);
    
    if plot_density
        figure;
        plot3k([ray_pos2(:,1), ray_pos2(:,3), nums], 'ColorBar', false);
        xlabel('x/mm');
        ylabel('y/mm');
        zlabel('Number of rays');
        set(gcf, 'PaperPosition', [0 0 14 12])
        set(gcf, 'Papersize', [14 12]);
        saveas(gcf, [thePath '/intensity_of_rays'], 'epsc')
    end
    
    % Repeat the starting positions
    f1 = @(k) repmat(ray_pos(k,:), round(nums(k)), 1);
    ray_pos = cell2mat(arrayfun(f1, (1:length(nums))', 'UniformOutput', false));
    n_rays = size(ray_pos, 1);
    
    ray_dir = repmat([1/sqrt(2), -1/sqrt(2), 0], n_rays, 1);
end

