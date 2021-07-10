% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
% RectangleInfo.m
% 
% Contains output data from a rectangular scan.
%
% PROPERTIES:
%  counters         - Contains the number of detected rays that had undergone
%                     1,2,3,4,etc. scatters for each pixel
%  cntrSum          - A matrix of the number of all detected rays for each pixel
%  counter_effusive - A matrix of the number of detected rays from the effuse
%                     beam
%  nx_pixels        - The number of pixels in the x direction
%  ny_pixels        - The number of pixels in the z direction
%  N_pixels         - The total number of pixels
%  rays_per_pixel   - The number of rays that were used per pixel
%  sample_surface   - A TriagSurface object of the sample used in the
%                     simulation
%  xrange           - The range of the scan, [min max]/mm, in the x direction
%  zrange           - The range of the scan, [min max]/mn, in the z direction
%  raster_movment   - The movment of the sample, in mm, between each pixel
%  time             - the time in seconds the simulation took
%  time_estimate    - the initial estimate of how long the simulation
%                     would take, in seconds
%
% METHODS:
%  TODO
classdef RectangleInfo < SimulationInfo

    properties %(SetAccess = immutable)
        n_detector;
        cntrSum;
        counter_effusive;
        num_killed;
        nx_pixels;
        nz_pixels;
        N_pixels;
        n_effuse;
        sample_surface;
        xrange;
        zrange;
        raster_movment_x;
        raster_movment_z;
        raster_pattern;
    end % End properties
    
    methods
        function obj = RectangleInfo(counters, num_killed, sample_surface, ...
                xrange, zrange, raster_movment_x, raster_movment_z, ...
                rays_per_pixel, n_effuse, time, t_estimate, cntr_effuse, ...
                n_detector, maxScatter, dist_to_sample, direct_beam, raster_pattern)
            if nargin ~= 17
                error('Wrong numer of input arguments');
            end
            obj = obj@SimulationInfo(time, t_estimate, direct_beam.init_angle, direct_beam, ...
                dist_to_sample, rays_per_pixel);
            obj.n_detector = n_detector;
            obj.num_killed = num_killed;
            obj.nx_pixels = raster_pattern.nx;
            obj.nz_pixels = raster_pattern.nz;
            obj.N_pixels = obj.nx_pixels*obj.nz_pixels;
            obj.sample_surface = sample_surface;
            obj.xrange = xrange;
            obj.zrange = zrange;
            obj.counters = {};
            for i_=1:obj.n_detector
                obj.counters{i_} = reshape(counters(:,i_,:,:), maxScatter, ...
                    obj.nz_pixels, obj.nx_pixels);
            end
            obj.raster_movment_x = raster_movment_x;
            obj.raster_movment_z = raster_movment_z;
            cntrSum = sum(counters,1);
            cntrSum2 = zeros(obj.n_detector, obj.nz_pixels, obj.nx_pixels);
            for k=1:obj.n_detector
                for i_=1:obj.nx_pixels
                    for j_=1:obj.nz_pixels
                        cntrSum2(k, j_, i_) = cntrSum(1, k, j_, i_);
                    end
                end
            end
            cntrSum = cntrSum2 + cntr_effuse;
            for i_=1:obj.n_detector
                obj.cntrSum{i_} = reshape(cntrSum(i_,:,:), obj.nz_pixels, ...
                    obj.nx_pixels);
            end
            for i_=1:obj.n_detector
                obj.counter_effusive{i_} = reshape(cntr_effuse(i_,:,:), obj.nz_pixels, ...
                    obj.nx_pixels);
            end
            obj.n_effuse = n_effuse;
            obj.raster_pattern = raster_pattern;
        end % End constructor
        
        function cntrSum2 = getSingle(obj, detector)
        % Gets the single scattering image matrix from the results object for
        % the specified detector.
        % 
        % Calling syntax:
        %  cntrSum = obj.getSingle(detector)
        %
        % INPUTS:
        %  detector - which detector to use, defaults to 1
        %
        % OUTPUT:
        %  cntrSum - raw image matrix of single scattering counts
            if nargin == 1
                detector = 1;
            end
            
            counters = obj.counters{detector}; %#ok<PROPLC>
            cntrSummed = counters(1,:,:); %#ok<PROPLC>
            cntrSum2 = zeros(obj.nz_pixels, obj.nx_pixels);
            for i_=1:obj.nx_pixels
                for j_=1:obj.nz_pixels
                    cntrSum2(j_,i_) = cntrSummed(1, j_, i_);
                end
            end
        end
        
        function cntrSum2 = getMultiple(obj, detector)
        % Gets the multiple scattering image matrix from the results object for
        % the specified detector.
        % 
        % Calling syntax:
        %  cntrSum = obj.getMultiple(detector)
        %
        % INPUTS:
        %  detector - which detector to use, defaults to 1
        %
        % OUTPUT:
        %  cntrSum - raw image matrix of multiple scattering counts
            if nargin == 1
                detector = 1;
            end
            
            counters2 = obj.counters{detector}(2:end, :, :);
            cntrSummed = sum(counters2, 1);
            cntrSum2 = zeros(obj.nz_pixels, obj.nx_pixels);
            for i_=1:obj.nx_pixels
                for j_=1:obj.nz_pixels
                    cntrSum2(j_,i_) = cntrSummed(1, j_, i_);
                end
            end
        end
        
        function addDetectorInfo(obj, aperture_c, aperture_axes)
            obj.aperture_c = aperture_c;
            obj.aperture_axes = aperture_axes;
        end
        
        function maxScatter = getMaxScatter(obj)
        % Gets the maximum number of scattering events that were allowed in the
        % simulation.
            maxScatter = size(obj.counters{1}, 1);
        end
        
        function cnts = countScattering(obj, n, detector)
        % Gives the nth scattering contribution (i.e. single scattering double
        % scattering etc.) or the range of scattering if n is provided as a 2
        % element vector. If the second limit is Inf then will count up to
        % infinity
            if nargin == 2
                detector = 1;
            end
            
            if length(n) == 1
                lims = n;
            elseif length(n) == 2
                if n(2) == Inf
                    lims = n(1):obj.getMaxScatter;
                else
                    lims = n(1):n(2);
                end
            else
                error('n must be of length 1 or 2');
            end
            
            cnts = reshape(sum(obj.counters{detector}(lims,:,:), 1), obj.nz_pixels, ...
                obj.nx_pixels);
        end
        
        function I = imageAll(obj, varargin)
        % Constructs an image based on the data in the RectangleInfo
        % object. Uses all contributions to the contrast, including the
        % effuse beam.
        %
        % INPUTS (all optional):
        %  scale        - Defines how to scale the gradient in the image,
        %                 can be: 'auto' - scaled with black being the
        %                           pixel with the lowest number of 
        %                           detection and white being the highest
        %                         'zeroed' - Black is made to be zero
        %                           detections with white still being the
        %                           pixel with the highest number
        %                         'removeZero; - black is the smallest
        %                           number of detection in a pixel that 
        %                           isn't zero
        %                         'manual' - The scale is specified manually
        %  specifyScale - The scale if it is being specified, [min max], if
        %                 'manual' is not used and image limits are then
        %                 this entry must be provided, but is ignored
        %  limX         - The x range of image to produce, If the limits
        %                 are not specifed they are automatically chosen
        %  limY         - The y range of image to produce, this corresponds
        %                 to the z direction in the simulation
        %  detector     - Which detector to use to produce the image
        %
        % OUTPUTS:
        %  I - A matrix of the grayscale data produced. The matrix is the
        %      size of the image and the value of each element is in [0,1] and
        %      gives the relative brightness of the associated pixel.
            
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'detector'
                        detector = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
            
            % Input checking
            if ~exist('scale', 'var')
                scale = 'auto';
            end
            if ~exist('specifyScale', 'var') && strcmp(scale, 'manual')
                error('Must specify a scale if you select manual scale.');
            end
            if ~exist('detector', 'var')
                detector = 1;
            end
            
            if exist('limX', 'var') && exist('limY', 'var')
                if strcmp('scale', 'manual')
                    I = obj.generalImage('im', obj.cntrSum{detector}, 'scale', ...
                        scale, 'specifyScale', specifyScale, 'limX', limX, 'limY', limY);
                else
                    I = obj.generalImage('im', obj.cntrSum{detector}, 'scale', ...
                        scale, 'limX', limX, 'limY', limY);
                end
            else
                if strcmp(scale, 'manual')
                    I = obj.generalImage('im', obj.cntrSum{detector}, 'scale', scale, 'specifyScale', specifyScale);
                else
                    I = obj.generalImage('im', obj.cntrSum{detector}, 'scale', scale);
                end
            end
        end % End image generating function
        
        function I = imageBigger(obj, scale)
        % TODO: update
        %
        % Artificially produce a larger image (2x) from the data. Can be useful
        % for small scans. The intermediate pixels are produced by
        % averaging the counts for the adjacent pixels.
        % 
        % One optional input:
        %  scale - string, 'auto' or 'zeroed'
            if nargin == 1
                scale = 'auto';
            end
            
            cntrSum2 = zeros(size(obj.cntrSum,1), 2*size(obj.cntrSum,2));
            
            for i_=1:size(obj.cntrSum,2)
                cntrSum2(:,2*i_) = obj.cntrSum(:,i_);
            end
            
            cntrSum3 = zeros(2*size(obj.cntrSum));
            
            for i_=1:size(obj.cntrSum,1)
                cntrSum3(2*i_,:) = cntrSum2(i_,:);
            end
            
            for i_=1:(0.5*size(cntrSum3, 2) - 1)
                ind = 2*i_;
                cntrSum3(:,ind+1) = (cntrSum3(:,ind) + cntrSum3(:,ind+2))/2;
            end
            
            for i_=1:(0.5*size(cntrSum3, 1) - 1)
                ind = 2*i_;
                cntrSum3(ind+1,:) = (cntrSum3(ind,:) + cntrSum3(ind+2,:))/2;
            end
            
            cntrSum3 = cntrSum3(2:end,2:end);
            
            if strcmp(scale, 'auto')
                I = mat2gray(cntrSum3);
            elseif strcmp(scale, 'zeroed')
                I = mat2gray(cntrSum3, [0 max(max(cntrSum3))]);
            else
                error('Specify a correct scaling mechanism');
            end
            
            I = I(1:end-1,1:end-1);
            RI = imref2d(size(I));
            RI.XWorldLimits = obj.xrange;
            RI.YWorldLimits = obj.zrange;
            figure
            imshow(I, RI);
            xlabel('-x/mm')
            ylabel('z/mm')
        end % End image enlarging function
        
        function [I, cntrSum2] = imageNoEffuse(obj, varargin)
        % Produces a 2d image from the data explicitly without the effuse
        % contribution. See imageAll.
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'detector'
                        detector = varargin{i_+1};
                    case 'plot'
                        make_plot = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
        
            % Input checking
            if ~exist('scale', 'var')
                scale = 'auto';
            end
            if ~exist('specifyScale', 'var') && strcmp(scale, 'manual')
                error('Must specify a scale if you select manual scale.');
            end
            if ~exist('detector', 'var')
                detector = 1;
            end
            if ~exist('make_plot', 'var')
                make_plot = false;
            end
            
            cntrSum2 = obj.cntrSum{detector} - obj.counter_effusive{detector};
            
            if exist('limX', 'var') && exist('limY', 'var')
                if strcmp('scale', 'manual')
                    I = obj.generalImage('im', cntrSum2, 'scale', scale, ...
                        'specifyScale', specifyScale, 'limX', limX, ...
                        'plot', make_plot, 'limY', limY);
                else
                    I = obj.generalImage('im', cntrSum2, 'scale', scale, ...
                        'limX', limX, 'limY', limY, 'plot', make_plot);
                end
            else
                I = obj.generalImage('im', cntrSum2, 'scale', scale, 'plot', ...
                    make_plot);
            end
        end
        
        function [I, cntrSum2] = imageSingle(obj, varargin)
        % Produces an image from the single scattering contribution only.
        % See imageAll.
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'detector'
                        detector = varargin{i_+1};
                    case 'plot'
                        make_plot = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
        
            % Input checking
            if ~exist('scale', 'var')
                scale = 'auto';
            end
            if ~exist('specifyScale', 'var') && strcmp(scale, 'manual')
                error('Must specify a scale if you select manual scale.');
            end
            if ~exist('detector', 'var')
                detector = 1;
            end
            if ~exist('make_plot', 'var')
                make_plot = true;
            end
            
            cntrSum2 = obj.getSingle(detector);
            
            if exist('limX', 'var') && exist('limY', 'var')
                if strcmp('scale', 'manual')
                    I = obj.generalImage('im', cntrSum2, 'scale', scale, ...
                        'specifyScale', specifyScale, 'limX', limX, ...
                        'plot', make_plot, 'limY', limY);
                else
                    I = obj.generalImage('im', cntrSum2, 'scale', ...
                        scale, 'limX', limX, 'limY', limY, 'plot', make_plot);
                end
            else
                I = obj.generalImage('im', cntrSum2, 'scale', scale, 'plot', ...
                    make_plot);
            end
        end
        
        function [I, cntrSum2] = imageMultiple(obj, varargin)
        % Produces a 2d image from the multiple scattering contribution
        % only. See imageAll.
        
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'detector'
                        detector = varargin{i_+1};
                    case 'plot'
                        make_plot = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
        
            % Input checking
            if ~exist('scale', 'var')
                scale = 'auto';
            end
            if ~exist('specifyScale', 'var') && strcmp(scale, 'manual')
                error('Must specify a scale if you select manual scale.');
            end
            if ~exist('detector', 'var')
                detector = 1;
            end
            if ~exist('make_plot', 'var')
                make_plot = true;
            end
            
            cntrSum2 = obj.getMultiple(detector);
            
            if exist('limX', 'var') && exist('limY', 'var')
                if strcmp('scale', 'manual')
                    I = obj.generalImage('im', cntrSum2, 'scale', scale, ...
                        'specifyScale', specifyScale, 'limX', limX, ...
                        'plot', make_plot, 'limY', limY);
                else
                    I = obj.generalImage('im', cntrSum2, 'scale', scale, ...
                        'limX', limX, 'limY', limY, 'plot', make_plot);
                end
            else
                I = obj.generalImage('im', cntrSum2, 'scale', scale, 'plot', ...
                    make_plot);
            end
        end
        
        function I = imageEffuse(obj, varargin)
        % Produces a 2d image from only the effusive contribtion. See
        % imageAll.
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'detector'
                        detector = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
            
            % Input checking
            if ~exist('scale', 'var')
                scale = 'auto';
            end
            if ~exist('specifyScale', 'var') && strcmp(scale, 'manual')
                error('Must specify a scale if you select manual scale.');
            end
            if ~exist('detector', 'var')
                detector = 1;
            end
            
            if exist('limX', 'var') && exist('limY', 'var')
                if strcmp('scale', 'manual')
                    I = obj.generalImage('im', obj.counter_effusive{detector}, 'scale', ...
                        scale, 'specifyScale', specifyScale, 'limX', limX, 'limY', limY);
                else
                    I = obj.generalImage('im', obj.counter_effusive{detector}, 'scale', ...
                        scale, 'limX', limX, 'limY', limY);
                end
            else
                I = obj.generalImage('im', obj.counter_effusive{detector}, 'scale', scale);
            end
        end
        
        function I = imageExtraEffuse(obj, varargin)
        % Produces an image with the effusive contribution enhanced.
        % INPUTS:
        %  factor - the factor to multiply the effuse contribution by
        %  scale  - how to scale the grayscale, see imageAll
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'detector'
                        detector = varargin{i_+1};
                    case 'factor'
                        factor = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
            
            if ~exist('scale', 'var')
                scale = 'auto';
            end
            if ~exist('specifyScale', 'var') && strcmp(scale, 'manual')
                error('Must specify a scale if you select manual scale.');
            end
            if ~exist('detector', 'var')
                detector = 1;
            end
            if ~exist('factor', 'var')
                factor = 1;
                warning('Using imageExtraEffuse without specifying how much extra');
            end
            cntrSum2 = obj.cntrSum{detector} + obj.counter_effusive{detector}*factor;
            
            if exist('limX', 'var') && exist('limY', 'var')
                if strcmp('scale', 'manual')
                    I = obj.generalImage('im', cntrSum2, 'scale', ...
                        scale, 'specifyScale', specifyScale, 'limX', limX, 'limY', limY);
                else
                    I = obj.generalImage('im', cntrSum2, 'scale', ...
                        scale, 'limX', limX, 'limY', limY);
                end
            else
                I = obj.generalImage('im', cntrSum2, 'scale', scale);
            end
        end
        
        function imageEyeView(obj, bar_location, Scalebar_length, dataPath)
            % TODO: parameters and default values etc.
            % Shifts on the scalebar from the edge of the image
            vert_shift = 0.07;
            hori_shift = 0.08;
            save_name = 'scale_var_all';

            % The size of the label font, pt
            font_size = 14;
            
            % Rotate the image
            I = rot90(mat2gray(obj.cntrSum), 3);
            x_lims = obj.zrange;
            y_lims = obj.xrange;
            x_res = range(x_lims)/obj.nz_pixels;
            y_res = range(y_lims)/obj.nx_pixels;
            
            
            % How many times larger are the vertical pixels than the horizontal pixels
            repeat_factor = y_res/x_res;
            if repeat_factor ~= 1
                if mod(repeat_factor, 1) ~= 0
                    non_int_scale = true;
                    warning('Non integer scaling between axis resolutions, ignoring.');
                else
                    non_int_scale = false;
                    tmp = [];
                    for i_=1:size(I, 1)
                        tmp = [tmp; repmat(I(i_,:), repeat_factor, 1)]; %#ok<AGROW>
                    end
                    I = tmp;
                end
            else
                non_int_scale = false;
            end

            figure;
            if non_int_scale
                y_lims = y_lims/repeat_factor;
                imagesc(x_lims, y_lims, I);
            else
                imagesc(x_lims, y_lims, I);
            end
            colormap(gray);
            axis equal tight
            hold on
            % Get the position of the scale bar
            switch bar_location
                case 'bottom left'
                    x_location = x_lims(2) - (1 - hori_shift)*abs(diff(x_lims)) + Scalebar_length;
                    y_location = y_lims(1) + (1 - 2*vert_shift)*abs(diff(y_lims));
                case 'bottom right'
                    x_location = x_lims(2) - hori_shift*abs(diff(x_lims));
                    y_location = y_lims(1) + (1 - 2*vert_shift)*abs(diff(y_lims));
                case 'top left'
                    x_location = x_lims(2) - (1 - hori_shift)*abs(diff(x_lims)) + Scalebar_length;
                    y_location = y_lims(1) + vert_shift*abs(diff(y_lims));
                case 'top right'
                    x_location = x_lims(2) - hori_shift*abs(diff(x_lims));
                    y_location = y_lims(1) + vert_shift*abs(diff(y_lims));
                otherwise
                    close;
                    error('Specify a correct location for the scale bar')
            end
            q_h = quiver(x_location, y_location, -Scalebar_length, 0, 'ShowArrowHead', ...
                        'off', 'AutoScale', 'off', 'LineWidth', 5, 'Color', 'y');
            bar_label = [num2str(Scalebar_length*1e3), '$\mu m$'];
            x_pos = x_location - Scalebar_length/2;
            y_pos = y_location + (font_size/20)*(vert_shift)*abs(diff(y_lims));
            t_h = text(x_pos, y_pos, bar_label, 'Color', 'y', 'FontSize', font_size, ...
                'HorizontalAlignment', 'center', 'Interpreter', 'latex');
            
            % Clear the axes
            axis off
            
            % Get the image size and the handles to the figure and axes
            im_size = size(I);
            fig_h = gcf;
            ax_h = gca;
            
            % Clear the greyspace
            set(ax_h, 'units', 'pixels') % set the axes units to pixels
            x = get(ax_h, 'position'); % get the position of the axes
            set(fig_h, 'units', 'pixels') % set the figure units to pixels
            y = get(fig_h, 'position'); % get the figure position
            set(fig_h, 'position', [y(1) y(2) (im_size(2)/im_size(1))*x(4) x(4)])% set the position of the figure to the length and width of the axes
            set(ax_h, 'units', 'normalized', 'position', [0 0 1 1]) % set the axes units to pixels

            % Set paper size for printing to pdf
            set(fig_h, 'Units', 'Inches');
            pos = get(fig_h, 'Position');
            set(fig_h, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)])

            % Save the figure to file

            % pdf
            print([dataPath '/' save_name '.eps'], '-depsc', '-r0')

            % fig
            savefig([dataPath '/' save_name '.fig'])
        end
        
        function contourImage(obj, n, limX, limY)
        % Produces a contour plot of the 2d data.
        %
        % INPUTS (all optional):
        %  n    - the number of contour levels to plot
        %  limX - the x limits of the plot
        %  limY - the y limits of the plot
            m = obj.cntrSum;
            m = flipud(m);
            xs = linspace(obj.xrange(1), obj.xrange(2), obj.nx_pixels);
            zs = linspace(obj.zrange(1), obj.zrange(2), obj.nz_pixels);
            [xs2, zs2] = meshgrid(xs, zs);
            figure
            if nargin == 1
                contour(xs2, zs2, m, 20);
            else
                contour(xs2, zs2, m, n);
            end
            axis('equal')
            if nargin > 3
                xlim(limX);
                ylim(limY);
            end
            xlabel('x/mm')
            ylabel('y/mm')
        end % End contour plot function
        
        function produceImages(obj, thePath)
        % Produces and saves a series of images from the simulation. Saves them
        % to the path provided. The colour scale takes black to be the lowest
        % number of counts and white to be the highest number of counts.
        %
        % Calling syntax:
        %  obj.produceImages(thePath)
        %
        % INPUT:
        %  thePath - data directory to save the images to
        
            for i_=1:obj.n_detector
                I = obj.imageSingle('detector', i_);
                title(['Detector ' num2str(i_)]);pause(0.1);
                imwrite(I, [thePath '/single' num2str(i_) '.png']);
                I = obj.imageMultiple('detector', i_);
                title(['Detector ' num2str(i_)]);pause(0.1);
                imwrite(I, [thePath '/multiple' num2str(i_) '.png']);
                if obj.n_effuse == 0
                    I = obj.imageAll('detector', i_);
                    title(['Detector ' num2str(i_)]);pause(0.1);
                    imwrite(I, [thePath '/noEffuse' num2str(i_) '.png']);
                else
                    I = obj.imageAll('detector', i_);
                    title(['Detector ' num2str(i_)]);pause(0.1);
                    imwrite(I, [thePath '/all' num2str(i_) '.png']);
                    I = obj.imageEffuse('detector', i_);
                    title(['Detector ' num2str(i_)]);pause(0.1);
                    imwrite(I, [thePath '/effuse' num2str(i_) '.png']);
                    I = obj.imageNoEffuse('detector', i_);
                    title(['Detector ' num2str(i_)]);pause(0.1);
                    imwrite(I, [thePath '/noEffuse' num2str(i_) '.png']);
                end
            end

            %obj.contourImage(30);
            %isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
            %if ~isOctave
            %    saveas(gcf, [thePath '/contourAll.eps'], 'epsc');
            %end
        end % End image creation and saving function.
        
        function ProduceImagesScaled(obj, thePath)
        % Produces and saves a series of images from the simulation, uses a
        % constant greyscale across all images
            
        end
        
        function counts = saveCounts(obj, fname, contribution, detector)
        % Saves a matrix of the counts of one of the contributions to a comma
        % delimated text file of the given name.
            if nargin == 3
                detector = 1;
            end
            counter = obj.counters{detector};
            switch contribution
                case 'single'
                    counts = squeeze(counter(1,:,:));
                case 'multiple'
                    counts = squeeze(sum(counter(2:end,:,:)));
                case 'effuse'
                    counts = obj.effuse_counters{detector};
            end
            dlmwrite(fname, counts, ',');
        end
        
        function saveText(obj, fname)
        % Saves a selection of theta to a text file for external plotting
        % and analysis. The text file is comma delinated and has the name
        % provided including the extension.
            fid = fopen(fname, 'w');
            
            fprintf(fid, '%s\n', 'x,y,Single,Multiple,Effuse');
            
            FORMAT = '%2.8f, ';
            
            % Create x,y,signal matrix for all the data
            
            % Wrtie the data to the text file
            
        end
        
        function [im, param, beam_param] = formatOutput(obj, dataPath)
        % Formats the results of a simulation into a more useful/simple format.
        % Useful for then using with Photo-Stereo reconstruction or similar. Ouputs
        % the data in structures rather than objects so that the class file is not
        % needed.
        % 
        % Calling syntax:
        %  [im, param, beam_param] = formatOutput(simData, dataPath)
        % 
        % INPUTS:
        %  simData  - The simulation data file, needs to be a rectangular scan
        %             result
        %  dataPath - Relative path to save the formatted output to. Output is put
        %             in a subdirectory in this path
        %
        % OUTPUTS:
        %  im         - Image results
        %  param      - General parameters of the simulation
        %  beam_param - Parameters for the set up of the beam
        
            % Core information on the image produced
            if obj.n_detector == 1
                [~, im.single{1}] = obj.imageSingle('plot', false);
                [~, im.multiple{1}] = obj.imageMultiple('plot', false);
            else
                for i_=1:obj.n_detector
                    [~, im.single{i_}] = obj.imageSingle('detector', i_, 'plot', ...
                        false);
                    [~, im.multiple{i_}] = obj.imageMultiple('detector', i_, 'plot', ...
                        false);
                end
            end
            im.raster_movement_x = obj.raster_movment_x;
            im.raster_movement_y = obj.raster_movment_z;

            % Main simulation parameters
            % Do not have the detector parameteres if an stl model of the
            % pinhole plate was used.
            if ~isnan(obj.aperture_c)
                if obj.n_detector == 1
                    param.detector_position{1} = obj.aperture_c;
                    param.detector_axes{1} = obj.aperture_axes;
                    x = -obj.aperture_c(1);
                    y = obj.aperture_c(2);
                    z = obj.dist_to_sample;
                    param.detector_vector{1} = [x, y, z]/norm([x, y, z]);
                else
                    for i_=1:obj.n_detector
                        inds = (2*i_:(2*i_+1)) - 1;
                        param.detector_position{i_} = obj.aperture_c(inds);
                        param.detector_axes{i_} = obj.aperture_axes(inds);
                        x = -obj.aperture_c(inds(1));
                        y = obj.aperture_c(inds(2));
                        z = obj.dist_to_sample;
                        param.detector_vector{i_} = [x, y, z]/norm([x, y, z]);
                    end
                end
            end
            param.z_sample_to_detector = obj.dist_to_sample;
            param.rays_per_pixel = obj.rays_per_pixel;

            % Parameters of the incidence beam
            beam_param.init_angle = obj.init_angle;
            beam_param.source_size = obj.beam_param.theta_max;
            beam_param.source_model = obj.beam_param.source_model;
            beam_param.pinhole_r = obj.beam_param.pinhole_r;
            beam_param.pinhole_c = obj.beam_param.pinhole_c;

            % Save the formatted data to a subdirectory 
            if ~exist([dataPath '/formatted'], 'dir')
                mkdir([dataPath '/formatted'])
            end
            dfile = [dataPath '/formatted/reconstructionSimulation.mat'];
            disp(['Saving to ' dfile]);
            save(dfile, 'im', 'param', 'beam_param')
        end        
    end % End public methods
    
    methods (Access = private)
        function I = generalImage(obj, varargin)
        % A general image constructing function, not to be called directly.
        % This function is called by all the other image generating
        % functions that plot just one image.
            for i_=1:2:length(varargin)
                switch varargin{i_}
                    case 'im'
                        cntrSum2 = varargin{i_+1};
                    case 'scale'
                        scale = varargin{i_+1};
                    case 'specifyScale'
                        specifyScale = varargin{i_+1};
                    case 'limX'
                        limX = varargin{i_+1};
                    case 'limY'
                        limY = varargin{i_+1};
                    case 'plot'
                        make_plot = varargin{i_+1};
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
            end
            
            % Default inputs
            if ~exist('make_plot', 'var')
                make_plot = true;
            end
            
            % Generate the data
            if strcmp('auto', scale)
                I = mat2gray(cntrSum2);
            elseif strcmp('zeroed', scale)
                I = mat2gray(cntrSum2, [0 max(max(cntrSum2))]);
            elseif strcmp('removeZero', scale)
                ind = cntrSum2 == 0;
                cntr = cntrSum2(~ind);
                I = mat2gray(cntrSum2, [min(cntr) max(cntr)]);
            elseif strcmp('manual', scale)
                I = mat2gray(cntrSum2, specifyScale);
            else
                error('must specify a valid gray scaling method');
            end
            
            if make_plot
                % Produce image
                RI = imref2d(size(I));
                RI.XWorldLimits = obj.xrange;
                RI.YWorldLimits = obj.zrange;
                figure
                imshow(I, RI);
                axis tight
                if exist('limX', 'var')
                    xlim(limX);
                end
                if exist('limY', 'var')
                    ylim(limY);
                end
                xlabel('-x/mm')
                ylabel('z/mm')
            end
        end % End general imaging function
    end % End private methods
end
