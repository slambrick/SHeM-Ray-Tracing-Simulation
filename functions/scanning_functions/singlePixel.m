% singlePixel.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Performs a single pixel simulation.
%
% TODO: name value pair inputs
%
% INPUTS:
%
% OUTPUTS:
%  simulationData - a SinglePixelInfor object contating infromation about the
%                   simulation that has been run
function simulationData = singlePixel(sample_surface, direct_beam, ...
        maxScatter, pinhole_surface, effuse_beam, ...
        dist_to_sample, sphere, circle, thePath, save_to_text, ...
        plate_represent, thePlate, apertureAbstract)
    
    % Estimation of the time the simulation will take
    t_estimate = time_estimate('n_rays', direct_beam{1}, 'n_effuse', ...
        effuse_beam{1}, 'sample_surface', sample_surface, 'n_pixels', 1, ...
        'pinhole_model', pinhole_model);    
    tic
    ray
    % Generate the effuse beam rays
    [effuse_pos, effuse_dir] = makeEffuse('n_effusive', effuse_beam{1}, 'pinhole_c', ...
        effuse_beam{2}, 'pinhole_r', effuse_beam{3}, 'cosine_n', effuse_beam{4});
    effuse_rays = {effuse_pos, effuse_dir};
    
    % Generate the direct beam rays
    [ray_pos, ray_dir] = makeDirect('n_rays', direct_beam{1}, ...
        'pinhole_c', direct_beam{2}, 'pinhole_r', direct_beam{3}, ...
        'plot_starting', false, 'theta_max', ...
        direct_beam{4}, 'source_model', direct_beam{5}, 'init_angle', direct_beam{6}, ...
        'sigma', direct_beam{7});
    rays = {ray_pos, ray_dir};
    
    switch plate_represent
        case 'stl'
            % Trace effuse beam
            [effuseCnt, effuseKilled, effuseLeft, ~, ~, ~, ~] = ...
                traceRays('rays', rays, 'sample', sample, 'maxScatter', ...
                    maxScatter, 'plate', pinhole_surface, 'scan_pos', scan_pos, ...
                    'dist', dist_to_sample, 'sphere', sphere, 'circle', circle);
            
            % Trace direct beam
            [cntr, killed, left, final_pos, final_dir, numScattersRay, ~] = ...
                traceRays('rays', effuse_rays, 'sample', sample, 'maxScatter', ...
                    maxScatter, 'plate', pinhole_surface, 'scan_pos', scan_pos, ...
                    'dist', dist_to_sample, 'sphere', sphere, 'circle', circle);
        case 'abstract'
            % TODO
        case 'circle'
            % Trace direct beam
            [cnt, killed, left, final_pos, final_dir, numScattersRay, ~] = ...
                traceSimple('rays', rays, 'sample', sample, 'maxScatter', ...
                    maxScatter, 'plate', thePlate, 'scan_pos', scan_pos, ...
                    'dist', dist_to_sample, 'sphere', sphere, 'circle', circle);
            
            % Trace effuse beam
            [effuseCnt, effuseKilled, effuseLeft, ~, ~, ~, ~] = ...
                traceSimple('rays', effuse_rays, 'sample', sample, 'maxScatter', ...
                    maxScatter, 'plate', thePlate, 'scan_pos', scan_pos, ...
                    'dist', dist_to_sample, 'sphere', sphere, 'circle', circle);
    end
      
    t = toc;
    
    % Generate the results object
    simulationData = SinglePixelInfo(cntr, killed, left, numScattersRay, ...
        final_pos, final_dir, t, effuseCnt, effuseKilled, effuseLeft);
    
    % The actual time the simulation took
    fprintf('Actual time taken: %f s\n', t);
    hr = floor(t/(60^2));
    min = floor((t - hr*60*60)/60);
    if min == 60
        hr = hr + 1;
        min = 0;
    end
    fprintf('Which is: %i hr %2i mins\n\n', hr, min);
    
    if save_to_text
        % Save the final positions and directions to a text file.
        simulationData.saveDirectionPositions([thePath ...
            '/positionsDirections.csv']);
    end
end

