classdef roughSurf1D
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vertices % Vertex points
        normals  % Surface normals
        N        % Number of surface segments
    end
    
    methods
        % Creates a random surface according to Bently's method
        function obj = random_surf_gen(h_RMS, Dx, corr_len, N)
            if ~exist('N', 'var')
                N = 10001;
            end
            if mod(N, 2) == 0
                error('Number of points in the random surface must be odd')
            end
            
            lambd = 2*corr_len^(2/3);
            h = h_RMS*sqrt(tanh(Dx/lambd));
            Z = randn(N, 1);
            ms = linspacee(-round(N/2), round(N/2), N);
            hs = h*conv(Z, e, 'same');
            xs = ms*Dx;
            
            obj.vertices = [xs; hs];
            ns = [hs(1:end-1) - hs(2:end); abs(xs(2:end) - xs(1:end-1))];
            obj.normals = ns./sqrt(sum(ns.^2, 1));
            obj.N = length(obj.normals);
        end % End surface generation method
        
        % Reads a random surface from a text file
        function obj = load_2Dsurface(fname)
            opts = delimitedTextImportOptions("NumVariables", 2);

            % Specify range and delimiter
            opts.DataLines = [10, Inf];
            opts.Delimiter = ",";

            % Specify column names and types
            opts.VariableNames = ["x", "y"];
            opts.VariableTypes = ["double", "double"];
            opts.ExtraColumnsRule = "ignore";
            opts.EmptyLineRule = "read";

            % Import the data
            surfaceused = readtable(fname, opts);
            xs = surfaceused.x;
            hs = surfaceused.y;
            
            obj.vertices = [x; y];
            ns = [hs(1:end-1) - hs(2:end); abs(xs(2:end) - xs(1:end-1))];
            obj.normals = ns./sqrt(sum(ns.^2, 1));
            obj.N = length(obj.normals);
        end % End load surface method
        
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
        function [analytic_thetas, weights] = alalytic_random_scatter(obj, init_angle, make_plots)
            % The incidence direction
            init_dir = [-sind(init_angle); -cosd(init_angle)];

            % Final directions across the whole surface
            % TODO: vectorise this
            analytic_dirs = zeros(size(obj.normals));
            for i_=1:length(analytic_dirs)
                analytic_dirs(:,i_) = init_dir - 2*(dot(obj.normals(:,i_), init_dir))*obj.normals(:,i_);
            end

            % Tangent vectors
            tangents = [obj.xs(2:end) - obj.xs(1:end-1); obj.hs(2:end) - obj.hs(1:end-1)];
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
                    (sum(analytic_dirs(2,:) < 0))/(obj.N - 1));
            end
            weights(analytic_dirs(2,:) < 0) = 0;

            if sum(weights < 0) < obj.N/20
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
        end % End analytic scatter method
    end % End methods
end

