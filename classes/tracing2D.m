% Copyright (c) 2018, 2020, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the  
% GNU/GPL-3.0-or-later.The
classdef tracing2D
    
    properties
        % No properties, just grouping funcitons
    end
    
    methods
        % Plots a multiple scattering vs single scattering histogram
        function [f, propotion_multiple] = multiple_scattering_hist(num_scatters, ratios)
            ind = num_scatters == 1;
            single = sum(ind, 2);
            multiple = sum(~ind, 2);
            propotion_multiple = multiple./(single + multiple);

            f = figure;
            plot(ratios, propotion_multiple)
            ylim([0, 1])
            xlabel('RMS height/correlation length')
            ylabel('Proportion multiply scattered')
            title('Proporion of rays multiply scattered');
        end
        
        % Calculates the scattered intensity in direction theta for the provided
        % incidence angle and the provided width parameter. The result is normalised as
        % a PDF. Specify the angles in degrees.
        %
        % INPUTS:
        %  theta   - The angle to the surface normal to calculate the scattered
        %            sudo intensity for, deg
        %  theta_I - The incidence angle to the surface normal, deg
        %  sigma_P - The parameter standard deviation, specified the wdith of the
        %            distribution. It is the standard deviation of the Gaussian
        %            component of the distribution, *not* the standard deviation of the
        %            distribution.
        %
        % OUTPUT:
        %  intensity - The scattered intensity in direction theta
        function intensity = broad_specular_distribution(theta, theta_I, sigma_P)
            % The specular direction
            theta_S = -theta_I;

            % Normalisation
            N = integral(@(x) cosd(x).*exp(-(theta_S - x).^2/(2*sigma_P^2)), -90, 90);

            % Intensity
            intensity = (1./N).*cosd(theta).*exp(-(theta_S - theta).^2/(2*sigma_P^2));
        end
        
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
        
        % A function to calculate the mean and standard deviation of the broad specular
        % distribution given the parameters of the distribution.
        %
        % INPUTS:
        %  sigmaP - The 'parameter standard deviation' of the distribution. This is the
        %           standard deviation of the Gaussian component of the distribution.
        %  thetaS - The specular angle to the normal.
        %
        % OUPUTS:
        %  mu    - The mean of the distibution 
        %  sigma - The standard deviation of the distribution.
        function [mu, sigma] = sigma_of_sigmaP(sigmaP, thetaS)
            % Normalisation of the pdf
            N = integral(@(x) cos(x).*normpdf(x, thetaS, sigmaP), -pi/2, pi/2);

            % Mean of the distribution
            mu = integral(@(x) x.*cos(x).*normpdf(x, thetaS, sigmaP), -pi/2, pi/2)/N;

            % Standard deviation
            sigma = integral(@(x) (x - mu).^2.*cos(x).*normpdf(x, thetaS, sigmaP)/N, -pi/2, pi/2);
            sigma = sqrt(sigma);
        end
        
        % Scatters rays uniformally off a randomly rough surface specified by the ratio
        % between the rms height and the correlation length. Provide the number of rays
        % and surface elements.
        %     The surface to be scattered off can either be specified by its statistics
        % or by explicity giving the x and y positions of the profile - do not give
        % both, the function will error if you do. Specify inputs as name-value pairs:
        % ('name', value).
        %
        % INPUTS:
        %  ratio      - The ratio between the 
        %  Nelements  - The number of 1D surface elements to model.
        %  xs         - Alternativly specify the surface to be scattered off explicitly,
        %               specify the x positions of the surfcace
        %  hs         - Specuify the height positions of the surface
        %  n_rays     - The number of rays to model.
        %  init_angle - The incident angle specified in degrees.
        %  make_plots - Boolean. Should plots be made? Plots of an example region of
        %               surface, the height distribution function of the generated
        %               surface, and histograms of the two components of the final
        %               directions and the number of scattering events. Optional,
        %               defaults to false.
        %  scattering - string, the type of scattering to perform
        %
        % OUTPUTS:
        %  thetas       - A vector of the final angles to the normal of the rays after
        %                 they have been scattered off the random surface.
        %  num_scatters - The number of scattering events that the rays have undergone.
        function [thetas, num_scatters] = random_scatter(varargin)
            % Get inputs
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'ratio'
                        ratio = varargin{i_+1};
                    case 'Nelements'
                        Nelements = varargin{i_+1};
                    case 'xs'
                        xs = varargin{i_+1};
                    case 'hs'
                        hs = varargin{i_+1};
                    case 'n_rays'
                        n_rays = varargin{i_+1};
                    case 'init_angle'
                        init_angle = varargin{i_+1};
                    case 'make_plots'
                        make_plots = varargin{i_+1};
                    case 'scattering'
                        scattering = varargin{i_+1};
                    case 'scattering_parameters'
                        scattering_parameters = varargin{i_+1};
                end
            end

            % Check that inputs exist and default values for non-existant inputs
            if ~exist('n_rays', 'var')
                error('Specify a number of rays.')
            end
            if ~exist('init_angle', 'var')
                error('Specify an incident angle.')
            end
            if ~exist('make_plots', 'var')
                make_plots = false;
            end

            % If the type of scattering isn't broad specular then the scattering
            % parameter isn't used
            if ~strcmp(scattering, 'broad specular')
                scattering_parameters = [0]; %#ok<NBRAK>
            elseif ~exist('scattering_parameters', 'var')
                error(['If using the broad specular distribution as the intrinsic ' ... 
                    'distribution then you must specify the broadness parameter.'])
            end

            % Switch between the cases where we are provided a surface profile and
            % provided the statistics on the surface profile
            if exist('ratio', 'var') && exist('Nelements', 'var')
                if exist('xs', 'var') || exist('hs', 'var')
                    error(['Do not specify both a surface profile and statistic for' ...
                        'gernerating a surface profile. Do one or the other.'])
                end
                % Generate a random Gaussian surface with a Gaussian height distribution
                % function and an exponential correlation length.
                [hs, xs] = roughSurf1D.rsgene(Nelements+1, (Nelements+1)/100, ratio, 1);
            elseif exist('xs', 'var') && exist('hs', 'var')
                if exist('ratio', 'var') || exist('Nelements', 'var')
                    error(['Do not specify both a surface profile and statistic for' ...
                        'gernerating a surface profile. Do one or the other.'])
                end
            else
                error(['Provide either a surface profile or the statistics for a' ...
                    'profile to be generated.'])
            end

            % Make plots if we want
            if make_plots
                % Plot a section of surface with equal axes
                figure
                plot(xs, hs)
                xlim([-ratio*20, ratio*20])
                axis equal
                title('Example surface profile')

                % Plot a histogram of the heights overlayed with the Gaussian distribution it
                % was generated from
                figure
                histogram(hs, 'Normalization', 'pdf')
                hold on
                height = linspace(min(hs), max(hs), 1000);
                plot(height, normpdf(height, 0, ratio))
                hold off
                xlabel('Height of surface in units of the correlation length')
                legend('Generated surface', 'Ideal surface')
            end

            % The minimum and maxium nominal intersection points
            min_x = min(xs)/2;
            max_x = max(xs)/2;

            % The starting positions
            init_y = range(xs);
            init_x = linspace(min_x, max_x, n_rays);
            init_pos = [init_x; repmat(init_y, 1, n_rays)];
            init_pos = init_pos - [tand(init_angle)*range(xs); 0];

            % The starting directions
            init_dir = [sind(init_angle); -cosd(init_angle)];

            % Calculate the normals to each surface elements
            vertices = [xs; hs];
            normals = [hs(1:end-1) - hs(2:end); abs(xs(2:end) - xs(1:end-1))];
            normals = normals./sqrt(sum(normals.^2, 1));

            % Scatter the rays using the ray tracing mex program
            [final_dirs, ~, num_scatters] = scatter_rays({vertices, normals}, ...
                {init_pos, init_dir}, {scattering, scattering_parameters});

            % Calculate the angle to the surface normal for the final directions of the rays
            thetas = atand(final_dirs(1,:)./final_dirs(2,:));

            if make_plots       
                % Plot a histogram of the final x and y componenets of the directions
                figure
                histogram(final_dirs(1,:), 'Normalization', 'pdf')
                xlabel('x component of direction')

                figure
                histogram(final_dirs(2,:), 'Normalization', 'pdf')
                xlabel('y component of direction')

                % Plot a histogram of the number of scattering events
                figure
                histogram(num_scatters, 'BinWidth', 1, 'Normalization', 'pdf')
                xlabel('Number of scatters')
            end
        end
    end
end

