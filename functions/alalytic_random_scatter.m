% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the Sub-beam Ray Tracing simulation, subject to the  
% GNU/GPL-3.0-or-later.
%
% analytic_random_scatter.m
%
% Uses a (more) analytic model to calculate the scattering distribution from a
% random surface -- assumes only single scattering. The final angle from an
% element of surface is calculated from the incident angle and the angle of the
% surface element, the weighting of the result from the element is calucated
% from the dot product of the element with the incidence direction.
%     Warnings are produced if it appears that the surface appears too rought
% for this model to be applied, although it is possible that no warning is
% produced and the model is not applicabel. Where the incidence angle (angle to
% the surface normal) is large the limit of the model is approached quicker (for
% less rough surfaces).
%
% INPUTS:
%  xs         - The x positions of the surface points
%  hs         - The heights of the surface points
%  init_angle - The incidence angle in degrees
%
% OUTPUTS:
%  analytic_thetas - The final directions from all the surface elements, in
%                    degrees.
%  weights         - The weighting 
function [analytic_thetas, weights] = alalytic_random_scatter(xs, hs, ...
        init_angle, make_plots)
    % The incidence direction
    init_dir = [-sind(init_angle); -cosd(init_angle)];
    
    % The number of surface points
    N = length(xs);
    
    % Calculate the normals to the surface elements
    normals = [hs(1:end-1) - hs(2:end); abs(xs(2:end) - xs(1:end-1))];
    normals = normals./sqrt(sum(normals.^2, 1));
    
    % Final directions across the whole surface
    % TODO: vectorise this
    analytic_dirs = zeros(size(normals));
    for i_=1:length(analytic_dirs)
        analytic_dirs(:,i_) = init_dir - 2*(dot(normals(:,i_), init_dir))*normals(:,i_);
    end
    
    % Tangent vectors
    tangents = [xs(2:end) - xs(1:end-1); hs(2:end) - hs(1:end-1)];
    tangents = tangents./sqrt(sum(tangents.^2, 1));
    
    % Weights of the directions
    weights = zeros(1, length(analytic_dirs));
    for i_ = 1:length(weights)
        weights(i_) = dot(init_dir, tangents(:,i_));
    end
    
    % Directions into the surface are given 0 weight
    if sum(analytic_dirs(2,:) < 0) > 0
        warning('Surface possibly too rough for single scattering model')
        fprintf('%f proportion of rays go into the surface\n', ...
            (sum(analytic_dirs(2,:) < 0))/(N - 1));
    end
    weights(analytic_dirs(2,:) < 0) = 0;
    
    if sum(weights < 0) < N/20
        if sum(weights < 0) > 0
            warning('Surface possibly too rough for single scattering model')
        end
        
        % Calculate the angles to the normal
        analytic_thetas = atand(analytic_dirs(1,:)./analytic_dirs(2,:));
        
        if make_plots
            % Plot of the resulting histogram
            weighted_histogram(analytic_thetas, weights, 75);
            xlabel('\theta/\circ')
            title('Single scattering off the surface')
        end
    else
        % Calculate the angles to the normal
        analytic_thetas = atand(analytic_dirs(1,:)./-analytic_dirs(2,:));
        warning(['Surface too rough for single scattering model over 5% or' ...
            'the results cannot be calculated.'])
    end
end

