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
classdef RectangleInfo

    properties %(SetAccess = immutable)
        counters;
        n_detector;
        cntrSum;
        counter_effusive;
        num_killed;
        nx_pixels;
        nz_pixels;
        N_pixels;
        rays_per_pixel;
        sample_surface;
        xrange;
        zrange;
        raster_movment_x;
        raster_movment_z;
        time;
        time_estimate;
    end % End properties
    
    methods
        function obj = RectangleInfo(counters, num_killed, sample_surface, ...
                xrange, zrange, raster_movment_x, raster_movment_z, ...
                rays_per_pixel, time, t_estimate, cntr_effuse, n_detector)
            if nargin ~= 12
                error('Wrong numer of input arguments');
            else
                obj.n_detector = n_detector;
                obj.num_killed = num_killed;
                obj.nx_pixels = length(xrange(1):raster_movment_x:xrange(2));
                obj.nz_pixels = length(zrange(1):raster_movment_z:zrange(2));
                obj.N_pixels = obj.nx_pixels*obj.nz_pixels;
                obj.sample_surface = sample_surface;
                obj.xrange = xrange;
                obj.zrange = zrange;
                for i_=1:obj.n_detector
                    obj.counters{i_} = reshape(counters(:,i_,:,:), maxScatter, ...
                        obj.nz_pixels, obj.nx_pixels);
                end
                obj.raster_movment_x = raster_movement_x;
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
                obj.rays_per_pixel = rays_per_pixel;
                obj.time = time;
                obj.time_estimate = t_estimate;
            end
        end % End constructor
        
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
                I = obj.generalImage('im', obj.cntrSum{detector}, 'scale', scale);
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
        
        function I = imageNoEffuse(obj, scale, specifyScale, limX, limY)
        % Produces a 2d image from the data explicitly without the effuse
        % contribution. See imageAll.
            cntrSum2 = obj.cntrSum - obj.counter_effusive;
            
            if nargin == 1
                I = obj.generalImage(cntrSum2, 'auto');
            elseif nargin == 2
                I = obj.generalImage(cntrSum2, scale);
            elseif nargin == 3
                I = obj.generalImage(cntrSum2, scale, specifyScale);
            elseif nargin == 5
                I = obj.generalImage(cntrSum2, scale, specifyScale, limX, limY);
            else
                error('Wrong number of arguments');
            end
        end
        
        function I = imageSingle(obj, scale, specifyScale, limX, limY)
        % Produces an image from the single scattering contribution only.
        % See imageAll.
            cntrSummed = obj.counters(1,:,:);
            cntrSum2 = zeros(obj.nz_pixels, obj.nx_pixels);
            for i_=1:obj.nx_pixels
                for j_=1:obj.nz_pixels
                    cntrSum2(j_,i_) = cntrSummed(1, j_, i_);
                end
            end
            
            if nargin == 1
                I = obj.generalImage(cntrSum2, 'auto');
            elseif nargin == 2
                I = obj.generalImage(cntrSum2, scale);
            elseif nargin == 3
                I = obj.generalImage(cntrSum2, scale, specifyScale);
            elseif nargin == 5
                I = obj.generalImage(cntrSum2, scale, specifyScale, limX, limY);
            else
                error('Wrong number of arguments');
            end
        end
        
        function I = imageMultiple(obj, scale, specifyScale, limX, limY)
        % Produces a 2d image from the multiple scattering contribution
        % only. See imageAll.
            counters2 = obj.counters(2:end, :, :);
            cntrSummed = sum(counters2, 1);
            cntrSum2 = zeros(obj.nz_pixels, obj.nx_pixels);
            for i_=1:obj.nx_pixels
                for j_=1:obj.nz_pixels
                    cntrSum2(j_,i_) = cntrSummed(1, j_, i_);
                end
            end
            
            if nargin == 1
                I = obj.generalImage(cntrSum2, 'auto');
            elseif nargin == 2
                I = obj.generalImage(cntrSum2, scale);
            elseif nargin == 3
                I = obj.generalImage(cntrSum2, scale, specifyScale);
            elseif nargin == 5
                I = obj.generalImage(cntrSum2, scale, specifyScale, limX, limY);
            else
                error('Wrong number of arguments');
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
            I = obj.imageSingle;
            imwrite(I,  [thePath '/single.png']);
            I = obj.imageMultiple;
            imwrite(I,  [thePath '/multiple.png']);
            obj.contourImage(30);
            isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
            if ~isOctave
                saveas(gcf, [thePath '/contourAll.eps'], 'epsc');
            end
            
            if sum(sum(obj.counter_effusive)) == 0
                I = obj.imageAll;
                imwrite(I, [thePath '/noEffuse.png']);
            else
                I = obj.imageAll;
                imwrite(I, [thePath '/all.png'])
                I = obj.imageNoEffuse;
                imwrite(I, [thePath '/noEffuse.png']);
                I = obj.imageEffuse;
                imwrite(I, [thePath '/effuse.png']);
            end
        end % End image creation and saving function.
        
        function counts = saveCounts(obj, fname, contribution)
        % Saves a matrix of the counts of one of the contributions to a comma
        % delimated text file of the given name.
            switch contribution
                case 'single'
                    counts = squeeze(obj.counters(1,:,:));
                case 'multiple'
                    counts = squeeze(sum(obj.counters(2:end,:,:)));
                case 'effuse'
                    counts = obj.effuse_counters;
            end
            dlmwrite(fname, counts, ',');
        end
    end % End methods
    
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
                    otherwise
                        warning(['Unknown input #' num2str(i_) ' to imageAll.']);
                end
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
            
            % Produce image
            RI = imref2d(size(I));
            RI.XWorldLimits = obj.xrange;
            RI.YWorldLimits = obj.zrange;
            figure
            imshow(I);%, RI);
            axis tight
            if specifylimits
                xlim(limX);
                ylim(limY);
            end
            xlabel('-x/mm')
            ylabel('z/mm')
        end % End general imaging function
    end % End private methods
end
