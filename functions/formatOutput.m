% formatOutput.m
%
% Copyright (c) 2019, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
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
function [im, param, beam_param] = formatOutput(simData, dataPath)
    % Core information on the image produced
    for i_=1:simData.n_detector
        [~, im.single{i_}] = simData.imageSingle('detector', i_, 'plot', ...
            false);
        [~, im.multiple{i_}] = simData.imageMultiple('detector', i_, 'plot', ...
            false);
    end
    im.raster_movement_x = simData.raster_movment_x;
    im.raster_movement_y = simData.raster_movment_z;
    
    % Main simulation parameters
    for i_=1:simData.n_detector
        inds = (2*i_:(2*i_+1)) - 1;
        param.detector_position{i_} = simData.aperture_c(inds);
        param.detector_axes{i_} = simData.aperture_axes(inds);
        x = simData.aperture_c(inds(1));
        y = simData.aperture_c(inds(2));
        z = simData.dist_to_sample;
        param.detector_vector{i_} = [x, y, z]/norm([x, y, z]);
    end
    param.z_sample_to_detector = simData.dist_to_sample;
    param.rays_per_pixel = simData.rays_per_pixel;
    
    % Parameters of the incidence beam
    beam_param.init_angle = simData.init_angle;
    beam_param.source_size = simData.beam_param.theta_max;
    beam_param.source_model = simData.beam_param.source_model;
    beam_param.pinhole_r = simData.beam_param.pinhole_r;
    beam_param.pinhole_c = simData.beam_param.pinhole_c;
    
    % Save the formatted data to a subdirectory 
    if ~exist([dataPath '/formatted'], 'dir')
        mkdir([dataPath '/formatted'])
    end
    dfile = [dataPath '/formatted/reconstructionSimulation.mat'];
    disp(['Saving to ' dfile]);
    save(dfile, 'im', 'param', 'beam_param')
end

