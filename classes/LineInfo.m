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
classdef LineInfo < SimulationInfo

    properties (SetAccess = private)
        Direction;
        sample_positions;
        range;
        N_pixels;
        single_scattering;
        multiple_scattering;
        num_killed;
        counters_effuse_single;
        counters_effuse_multiple;
        killed_effuse;
        raster_movment;
    end % End properties

    methods
        function obj = LineInfo(Direction, range, working_distance, counters, ...
                num_killed, raster_movment, n_rays, time, ...
                time_estimate, counters_effuse_single, ...
                counters_effuse_multiple, killed_effuse, direct_beam)
            obj = obj@SimulationInfo(time, time_estimate, direct_beam.init_angle, direct_beam, ...
                working_distance, n_rays);
            obj.Direction = Direction;
            obj.sample_positions = range(1):raster_movment:range(2);
            obj.range = range;
            obj.N_pixels = length(obj.sample_positions);
            obj.counters = counters;
            obj.single_scattering = counters(1,:);
            obj.multiple_scattering = sum(counters(2:end,:));
            obj.num_killed = num_killed;
            obj.raster_movment = raster_movment;
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
                obj.counters_effuse_multiple + obj.counters_effuse_single;
        end % End total function

        function [fig, ax] = linePlotParts(obj)
            % Produces a few plots of the line scan.
            single = obj.single_scattering;
            multiple = obj.multiple_scattering;
            effuse_single = obj.counters_effuse_single;
            effuse_multiple = obj.counters_effuse_multiple;
            xs = obj.sample_positions;

            [fig, ax] = obj.linePlot();
            hold(ax, 'on')
            plot(ax, xs, single, 'LineWidth', 3)
            plot(ax, xs, multiple, 'LineWidth', 3)
            
            % legend depending on whether effuse exists
            if sum(effuse_single + effuse_multiple) > 0
                plot(ax, xs, effuse_single + effuse_multiple, 'LineWidth', 3)
                legend(ax, 'Total', 'Single', 'Multiple', 'Effuse');
            else
                legend(ax, 'Total', 'Single', 'Multiple')
            end
            hold(ax, 'off')
        end
        
        function [fig, ax] = linePlot(obj)
            total = obj.totalLine();
            xs = obj.sample_positions;

            fig = figure;
            ax = axes(fig);
            plot(ax, xs, total, 'LineWidth', 3)

            if obj.Direction == 'y'
                xlabel(ax, 'Distance from pinhole plate/mm')
            else
                xlabel(ax, [obj.Direction '/mm']);
            end
            ylabel(ax, 'Counts');

            % size and font for report
            xlim(ax, [obj.sample_positions(1), obj.sample_positions(end)])
            set(ax, 'FontSize', 30)
            set(fig, 'Position', [100 100 900 800])
        end
        
        function producePlots(obj, thePath)
            [f1, ~] = obj.linePlotParts();
    
            if nargin == 2
                saveas(f1, [thePath '/line_comparison_plot.eps'], 'epsc');
            end

            if nargin == 2
                saveas(gcf, [thePath '/line_plot.eps'], 'epsc');
            end
        end % End plotting function
        
        function [im, param, beam_param] = formatOutput(obj, dataPath)
            % Core information on the image produced
            im.single = obj.single_scattering;
            im.multiple = obj.multiple_scattering;
            im.effuse = obj.counters_effuse_multiple + obj.counters_effuse_single;
            im.raster_movement = obj.raster_movment;
            param.sample_positions = obj.sample_positions;
            param.axis = obj.Direction;
            
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
    end

end

