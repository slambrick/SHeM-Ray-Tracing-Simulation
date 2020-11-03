% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% SinglePixelInfo.m
%
% Contains the results froma line scan simulation:
%
% PROPERTIES:
%  number_of_scatters_pre_ray - The number of scattering events that each
%                               detected ray underwent.
%  final_positions    - a matrix of the positions where the detected rays
%                       hit the detector surface.
%  final_directions   - a matrix of the directions the detected rays were
%                       travelling when they hit the detector surface.
%  number_killed      - the number of rays that had to be forcibly stopped
%                       because they exceeded the maximum number of scatters.
%  number_detected    - the number of detected rays
%  number_left        - the number of rays that left the region o
%                       simulation, they did not scatter off either the 
%                       pinhole plate or sample, or get detected.
%  time               - the time in seconds the simulation took
%  effuse_detected    - the number of rays from the effuse beam that were
%                       detected
%  effuse_killed      - the number of rays from the effuse beam that were
%                       forcibly stopped
%  effuse_left        - the number of rays from the effuse beam that left the
%                       simulation volume
classdef SinglePixelInfo

    properties (SetAccess = immutable)
        number_of_scatters_per_ray;
        final_positions;
        final_directions;
        number_killed;
        number_detected;
        number_left;
        time;
        effuse_detected;
        effuse_killed;
        effuse_left;
    end
    
    methods
        function obj = SinglePixelInfo(cntr, killed, left, numScattersRay, ...
                final_pos, final_dir, time, effuse_det, effuse_killed, ...
                effuse_left)
        % Constructor, arguments are the properties of the object.
            obj.number_of_scatters_per_ray = numScattersRay;
            obj.final_positions = final_pos;
            obj.final_directions = final_dir;
            obj.number_killed = killed;
            obj.number_detected = cntr;
            obj.number_left = left;
            obj.time = time;
            obj.effuse_detected = effuse_det;
            obj.effuse_killed = effuse_killed;
            obj.effuse_left = effuse_left;
        end % End constructor.
        
        function saveDirectionPositions(obj, fname)
        % Saves the final positions and directions of rays, along with the 
        % number of sample scatters each ray undergoes to a text file. 
        % Deliminates with a comma and saves to the provided filename (including
        % extension).
            fid = fopen(fname, 'w');
            
            fprintf(fid, '%s\n', 'Direction x, Direction y, Direction z, Position x, Position y, Position z, Number of scatters');
            
            n = length(obj.number_of_scatters_per_ray);
            FORMAT = '%2.8f, ';
            for i_=1:n
                fprintf(fid, FORMAT, obj.final_directions(i_, 1));
                fprintf(fid, FORMAT, obj.final_directions(i_, 2));
                fprintf(fid, FORMAT, obj.final_directions(i_, 3));
                fprintf(fid, FORMAT, obj.final_positions(i_, 1));
                fprintf(fid, FORMAT, obj.final_positions(i_, 2));
                fprintf(fid, FORMAT, obj.final_positions(i_, 3));
                fprintf(fid, '%i\n', obj.number_of_scatters_per_ray(i_));
            end
            
            fclose(fid);
        end % End outputting function
        
        function [single, multiple, effuse] = single_multi_effuse(obj)
        % Gives the number of singly scattered (sample), multiply scattered
        % and effusive rays detected for the given simulation.
            effuse = obj.effuse_detected;
            
            ind = obj.number_of_scatters_per_ray == 1;
            single = size(obj.number_of_scatters_per_ray(ind), 1);
            multiple = size(obj.number_of_scatters_per_ray, 1) - single;
        end % End simple binning function
        
        function histRays = scatteringHistogram(obj, mmax)
        % Produce a histogram of the number of sample scatters that the rays
        % underwent. Goes from 1 upto the value of mmax.
            histRays = binMyWayMex(obj.number_of_scatters_per_ray, mmax);
            figure
            xs = 1:mmax;
            bar(xs, histRays);
            xlabel('Number of sample scattering events');
            ylabel('Number of rays detected');
        end
    end % End methods
end

