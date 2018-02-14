% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

classdef RectangleInfo
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

    properties (SetAccess = immutable)
        counters;
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
        raster_movment;
        time;
        time_estimate;
    end % End properties
    
    methods
        function obj = RectangleInfo(counters, num_killed, sample_surface, ...
                                     xrange, zrange, raster_movment, ...
                                     rays_per_pixel, time, t_estimate, cntr_effuse)
            if nargin ~= 10
                error('Wrong numer of input arguments');
            else
                obj.counters = counters;
                obj.num_killed = num_killed;
                obj.nx_pixels = length(xrange(1):raster_movment:xrange(2));
                obj.nz_pixels = length(zrange(1):raster_movment:zrange(2));
                obj.N_pixels = obj.nx_pixels*obj.nz_pixels;
                obj.sample_surface = sample_surface;
                obj.xrange = xrange;
                obj.zrange = zrange;
                obj.raster_movment = raster_movment;
                cntrSum = sum(counters,1);
                cntrSum2 = zeros(obj.nz_pixels, obj.nx_pixels);
                for i_=1:obj.nx_pixels
                    for j_=1:obj.nz_pixels
                        cntrSum2(j_,i_) = cntrSum(1, j_, i_);
                    end
                end
                obj.cntrSum = cntrSum2 + cntr_effuse;
                obj.counter_effusive = cntr_effuse;
                obj.rays_per_pixel = rays_per_pixel;
                obj.time = time;
                obj.time_estimate = t_estimate;
            end
        end % End constructor
        
        function I = imageAll(obj, scale, specifyScale, limX, limY)
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
        %
        % OUTPUTS:
        %  I - A matrix of the grayscale data produced. The matrix is the
        %      size of the image and the value of each element is in [0,1] and
        %      gives the relative brightness of the associated pixel.
            
            if nargin == 1
                I = obj.generalImage(obj.cntrSum, 'auto');
            elseif nargin == 2
                I = obj.generalImage(obj.cntrSum, scale);
            elseif nargin == 3
                I = obj.generalImage(obj.cntrSum, scale, specifyScale);
            elseif nargin == 5
                I = obj.generalImage(obj.cntrSum, scale, specifyScale, limX, limY);
            else
                error('Wrong number of arguments');
            end
        end % End image generating function
        
        function I = imageBigger(obj, scale)
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
        
        function I = imageEffuse(obj, scale, specifyScale, limX, limY)
        % Produces a 2d image from only the effusive contribtion. See
        % imageAll.
            if nargin == 1
                I = obj.generalImage(obj.counter_effusive, 'auto');
            elseif nargin == 2
                I = obj.generalImage(obj.counter_effusive, scale);
            elseif nargin == 3
                I = obj.generalImage(obj.counter_effusive, scale, specifyScale);
            elseif nargin == 5
                I = obj.generalImage(obj.counter_effusive, scale, ...
                                     specifyScale, limX, limY);
            else
                error('Wrong number of arguments')
            end
        end
        
        function I = imageExtraEffuse(obj, factor, scale, specifyScale, limX, limY)
        % Produces an image with the effusive contribution enhanced.
        % INPUTS:
        %  factor - the factor to multiply the effuse contribution by
        %  scale  - how to scale the grayscale, see imageAll
            cntrSum2 = obj.cntrSum + obj.counter_effusive*factor;
            
            if nargin == 2
                I = obj.generalImage(cntrSum2, 'auto');
            elseif nargin == 3
                I = obj.generalImage(cntrSum2, scale);
            elseif nargin == 5
                I = obj.generalImage(obj.counter_effusive, scale, ...
                                     specifyScale, limX, limY);
            else
                error('Wrong number of arguments')
            end
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
            obj.imageSingle;
            saveas(gcf, [thePath '/single.png'], 'png');
            obj.imageMultiple;
            saveas(gcf, [thePath '/multiple.png'], 'png');
            obj.contourImage(30);
            saveas(gcf, [thePath '/contourAll.eps'], 'epsc');
            
            if sum(sum(obj.counter_effusive)) == 0
                obj.imageAll;
                saveas(gcf, [thePath '/noEffuse.png'], 'png');
            else
                obj.imageAll;
                saveas(gcf, [thePath '/all.png'], 'png');
                obj.imageNoEffuse;
                saveas(gcf, [thePath '/noEffuse.png'], 'png');
                obj.imageEffuse;
                saveas(gcf, [thePath '/effuse.png'], 'png');
            end
        end % End image creation and saving function.
    end % End methods
    
    methods (Access = private)
        function I = generalImage(obj, cntrSum2, scale, specifyScale, limX, limY)
        % A general image constructing function, not to be called directly.
        % This function is called by all the other image generating
        % functions that plot just one image.
            if nargin == 3
                specifylimits = false;
            elseif nargin == 4
                specifylimits = false;
            else
                specifylimits = true;
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
            imshow(I, RI);
            if specifylimits
                xlim(limX);
                ylim(limY);
            end
            xlabel('-x/mm')
            ylabel('z/mm')
        end % End general imaging function
    end % End private methods
end
