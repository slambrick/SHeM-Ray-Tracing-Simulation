% LineInfo.m
%
% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Contains the results from a line scan simulation.
%
% PROPERTIES:
%  Direction           - Was the line scan performed in the x, y or z direction
%  sample_positions    - The positions, in wichever direction the scan is
%                        taken, that the sample is rastered through
%  range               - The range of scane, [min max], in mm
%  N_pixels            - The number of pixels
%  counters            - Contains the number of detected rays that had undergone
%                        1,2,3,4,etc. scatters for each pixel
%  single_scattering   - A vector of the number of rays from the direct beam
%                        that were detected that scattered once off the sample
%  multiple_scattering - A vector of the number of rays from the direct beam
%                        that were detected that scattered more than once off 
%                        the sample
%  num_killed          - A vector of the number of rays that had to be 
%                        forcibly stopped because they exceeded the maximum 
%                        allowed number of scatters for each pixel
%  counters_effuse_single   - A vector of the number of detected rays from the
%                             effuse beam that scattered once off the sample
%  counters_effuse_multiple - A vector of the number of detected rays form teh
%                             effuse beam that scattered more than once off the 
%                             sample
%  killed_effuse       - A vector of the number of rays that were forcibly
%                        stopped from the effuse beam
%  raster_movment      - The movment of the sample, in mm, between each pixel
%  rays_per_pixel      - The number of rays that were used per pixel
%  time                - the time in seconds the simulation took
%  time_estimate       - the initial estimate of how long the simulation
%                        would take, in seconds
classdef LineInfo

    properties (SetAccess = private)
        Direction;
        sample_positions;
        range;
        N_pixels;
        counters;
        single_scattering;
        multiple_scattering;
        num_killed;
        counters_effuse_single;
        counters_effuse_multiple;
        killed_effuse;
        raster_movment;
        rays_per_pixel;
        time;
        time_estimate;
    end % End properties
    
    methods
        function obj = LineInfo(Direction, range, counters, ...
                                num_killed, raster_movment, n_rays, time, ...
                                time_estimate, counters_effuse_single, ...
                                counters_effuse_multiple, killed_effuse)
            obj.Direction = Direction;
            obj.sample_positions = range(1):raster_movment:range(2);
            obj.range = range;
            obj.N_pixels = length(obj.sample_positions);
            obj.counters = counters;
            obj.single_scattering = counters(1,:);
            obj.multiple_scattering = sum(counters(2:end,:));
            obj.num_killed = num_killed;
            obj.raster_movment = raster_movment;
            obj.rays_per_pixel = n_rays;
            obj.time = time;
            obj.time_estimate = time_estimate;
            obj.counters_effuse_single = counters_effuse_single;
            obj.counters_effuse_multiple = counters_effuse_multiple;
            obj.killed_effuse = killed_effuse;
        end % End constructor
        
        function saveText(obj, fname)
        % Saves a selection of the data to a text file for external plotting and
        % analysis. The text file is comma delinated and has the name provided,
        % including the extension.
            fid = fopen(fname, 'w');
            
            fprintf(fid, '%s\n', ['Positions,Single,Multiple,EffuseSingle,' ...
                'EffuseMultiple']);
            
            FORMAT = '%2.8f, ';
            
            for i_=1:length(obj.sample_positions)
                fprintf(fid, FORMAT, obj.sample_positions(i_));
                fprintf(fid, FORMAT, obj.single_scattering(i_));
                fprintf(fid, FORMAT, obj.multiple_scattering(i_));
                fprintf(fid, FORMAT, obj.counters_effuse_single(i_));
                fprintf(fid, '%2.8f\n', obj.counters_effuse_multiple(i_));
            end
            
            fclose(fid);
        end % End text file function
        
        function tot = totalLine(obj)
        % Gives the total contribution to the line scan.
            tot = obj.single_scattering + obj.multiple_scattering + ...
                obj.counters_effuse_multiple' + obj.counters_effuse_single';
        end % End total function
        
        function producePlots(obj, thePath)
        % Produces a few plots of the line scan.
            single = obj.single_scattering;
            multiple = obj.multiple_scattering;
            effuse_single = obj.counters_effuse_single';
            effuse_multiple = obj.counters_effuse_multiple';
            total = obj.totalLine;
            xs = obj.sample_positions;
            
            if obj.Direction == 'y'
                xs = xs + 2.121;
            end
            
            figure
            plot(xs, single)
            hold on
            plot(xs, multiple)
            plot(xs, total)
            if sum(effuse_single + effuse_multiple) > 0
                plot(xs, effuse_single + effuse_multiple)
                legend('Single', 'Multiple', 'Total', 'Effuse');
            else
                legend('Single', 'Multiple', 'Total')
            end
            
            if obj.Direction == 'y'
                xlabel('Distance from pinhole plate/mm')
            else
                xlabel([obj.Direction '/mm']);
            end
            ylabel('Number of rays detected');
            grid on
            if nargin == 2
                saveas(gcf, [thePath '/line_comparison_plot.eps'], 'epsc');
            end
            
            figure
            plot(xs, total)
            if obj.Direction == 'y'
                xlabel('Distance from pinhole plate/mm')
            else
                xlabel([obj.Direction '/mm']);
            end
            ylabel('Number of rays detected');
            grid on
            if nargin == 2
                saveas(gcf, [thePath '/line_plot.eps'], 'epsc');
            end
        end % End plotting function
    end
    
end

