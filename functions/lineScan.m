% lineScan.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates a 1d simulation of the sample.
%
% TODO: name value pair inputs
%
% Calling syntax:
%
% INPUTS:
%
% OUTPUTS:
%  line_scan_info - A LineInfo object containing the information about the
%                   simulation
function line_scan_info = lineScan(sample_surface, scan_range, direct_beam, ...
        raster_movement, maxScatter, Direction, pinhole_surface, effuse_beam, ...
        dist_to_sample, sphere, thePath, save_text, pinhole_model, thePlate, ...
        apertureAbstract, ray_model)
    
    % Dependent variables
    textFname = 'data_for_plotting.csv';
    
    % Sample positions
    sample_xs = scan_range(1):raster_movement:scan_range(2);
    n_pixels = length(sample_xs);
    
    % Create variables for output data
    counters = zeros(maxScatter, n_pixels);
    num_killed = zeros(n_pixels, 1);
    cntr_effuse_single = zeros(n_pixels, 1);
    counter_effuse_multiple = zeros(n_pixels, 1);
    killed_effuse = zeros(n_pixels, 1);
    
    % Estimate of the time for the simulation
    % TODO: change to estimate the time for the simple models
    t_estimate = time_estimate('n_rays', direct_beam.n, 'n_effuse', ...
        effuse_beam.n, 'sample_surface', sample_surface, 'n_pixels', n_pixels, ...
        'pinhole_model', pinhole_model);
    
    tic
    
    % Are we running in Matlab or Octave
    isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    
    % Starts the parallel pool if one does not already exist.
    if ~isOctave
        if isempty(gcp('nocreate'))
            parpool 
        end
    end
    
    
    % Generates a graphical progress bar if we are using the MATLAB GUI.
    if ~isOctave
        progressBar = feature('ShowFigureWindows');
    else
        progressBar = true;
    end
    if progressBar && ~isOctave
        ppm = ParforProgMon('Simulation progress: ', n_pixels);
        h = 0;
    elseif isOctave
        h = waitbar(0, 'Simulation progress: ');
    else 
        % If the variable ppm is undefined then the parfor loop will
        % throw errors.
        ppm = 0;
        h = 0;
    end
    
    % Makes the parfor loop stop complaining.
    plate_represent = pinhole_model;
    
    % TODO: make this parallel in Octave
    parfor i_=1:n_pixels
        scan_pos = sample_xs(i_);
        
        % Put the sample into the right place for this iteration
        switch Direction
            case 'x'
                this_surface = copy(sample_surface);
                this_surface.moveBy([scan_pos 0 0]);
                scan_pos_x = scan_pos;
                scan_pos_z = 0;
            case 'y'
                this_surface = copy(sample_surface);
                this_surface.moveBy([scan_pos -scan_pos 0]);
                scan_pos_x = scan_pos;
                scan_pos_z = 0;
            case 'z'
                this_surface = copy(sample_surface);
                this_surface.moveBy([0 0 scan_pos]);
                scan_pos_x = 0;
                scan_pos_z = scan_pos;
            otherwise
                error('Specify a correct direction for the line scan')
        end
        scan_pos2 = [scan_pos_x, scan_pos_z];
        
        % Direct beam
        [~, killed, numScattersRay] = switch_plate('plate_represent', ...
            plate_represent, 'sample', this_surface, 'maxScatter', maxScatter, ...
            'pinhole_surface', pinhole_surface, 'thePlate', thePlate, ...
            'dist', dist_to_sample, 'sphere', sphere, 'ray_model', ...
            ray_model, 'which_beam', direct_beam.source_model, 'beam', direct_beam);
        
        % Effuse beam
        [~, effuseKilled, numScattersEffuse] = switch_plate('plate_represent', ...
            plate_represent, 'sample', this_surface, 'maxScatter', maxScatter, ...
            'pinhole_surface', pinhole_surface, 'thePlate', thePlate, ...
            'dist', dist_to_sample, 'sphere', sphere, 'ray_model', ...
            ray_model, 'which_beam', 'Effuse', 'beam', effuse_beam);
        
        % Update the progress bar if we are working in the MATLAB GUI.
        if progressBar && ~isOctave
            ppm.increment();
        elseif isOctave
            waitbar(i_/n_pixels, h);
        end
        
        % Save the data for this iteration
        cntr_effuse_single(i_) = numScattersEffuse(1);
        counter_effuse_multiple(i_) = sum(numScattersEffuse(2:end));
        killed_effuse(i_) = effuseKilled;
        counters(:,i_) = numScattersRay;
        num_killed(i_) = killed;
        
        % Delete the surface object for this iteration
        if ~isOctave
            delete(this_surface);
        end
    end
    
    % Close the parallel pool
    current_pool = gcp('nocreate');
    delete(current_pool);
    
    t = toc;
    
    % The actual time the simulation took
    fprintf('Actual time taken: %f s\n', t);
    hr = floor(t/(60^2));
    min = floor((t - hr*60*60)/60);
    if min == 60
        hr = hr + 1;
        min = 0;
    end
    fprintf('Which is: %i hr %2i mins\n\n', hr, min);
    
    % Generate an object to store the information
    line_scan_info = LineInfo(Direction, scan_range, counters, num_killed, ...
        raster_movement, direct_beam.n, t, t_estimate, cntr_effuse_single, ...
        counter_effuse_multiple, killed_effuse);
    
    % Save desired outputs to text for plotting in other software
    if save_text
        line_scan_info.saveText([thePath '/' textFname]);
    end
    
    if progressBar
        line_scan_info.producePlots(thePath);
    end
end
