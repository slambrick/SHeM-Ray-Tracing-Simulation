% lineScan.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.
%
% Generates a 1d simulation of the sample.
%
% Calling syntax:
%
% INPUTS:
%
% OUTPUTS:
%  line_scan_info - A LineInfo object containing the information about the
%                   simulation
function line_scan_info = lineScan(varargin)
    for i_=1:2:length(varargin)
        switch varargin{i_}
            case 'sample_surface'
                sample_surface = varargin{i_+1};
            case 'scan_inputs'
                scan_inputs = varargin{i_+1};
            case 'direct_beam'
                direct_beam = varargin{i_+1};
            case 'max_scatter'
                maxScatter = varargin{i_+1};
            case 'pinhole_surface'
                pinhole_surface = varargin{i_+1};
            case 'effuse_beam'
                effuse_beam = varargin{i_+1};
            case 'dist_to_sample'
                dist_to_sample = varargin{i_+1};
            case 'sphere'
                sphere = varargin{i_+1};
            case 'circle'
                circle = varargin{i_+1};
            case 'thePath'
                thePath = varargin{i_+1};
            case 'pinhole_model'
                pinhole_model = varargin{i_+1};
            case 'thePlate'
                thePlate = varargin{i_+1};
            case 'aperture_abstract'
                aperture_abstract = varargin{i_+1};
            case 'ray_model'
                ray_model = varargin{i_+1};
            case 'n_detector'
                n_detector = varargin{i_+1};
            case 'make_plots'
                make_plots = varargin{i_+1};
            case 'init_angle'
                init_angle = varargin{i_+1};
            otherwise
                error(['input ' num2str(i_) ' not recognised:']);
        end
    end
    
    % Sample positions
    sample_xs = scan_inputs.range1D(1):scan_inputs.raster_movment1D:scan_inputs.range1D(2);
    n_pixels = length(sample_xs);

    % Create variables for output data
    counters = zeros(maxScatter, n_detector, n_pixels);
    num_killed = zeros(n_pixels, 1);
    cntr_effuse_single = zeros(n_detector, n_pixels);
    counter_effuse_multiple = zeros(n_detector, n_pixels);
    killed_effuse = zeros(n_detector, n_pixels);

    % Estimate of the time for the simulation
    % TODO: change to estimate the time for the simple models
    %t_estimate = time_estimate('n_rays', direct_beam.n, 'n_effuse', ...
    %    effuse_beam.n, 'sample_surface', sample_surface, 'n_pixels', n_pixels, ...
    %    'pinhole_model', pinhole_model);
    t_estimate = 0;
    
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
    progressBar = false

    ppm = 0;
    if progressBar
        if ~isOctave
            ppm = ParforProgressbar(n_pixels, 'showWorkerProgress', true);
            h = '';
        else
            h = waitbar(0, 'Simulation progress: ');
        end
    else
        h = '';
    end

    % Makes the parfor loop stop complaining.
    plate_represent = pinhole_model;

    % TODO: make this parallel in Octave
    parfor i_=1:n_pixels
        scan_pos = sample_xs(i_);
        this_surface = copy(sample_surface);

        % Put the sample into the right place for this iteration
        switch scan_inputs.direction_1D
            case 'x'
                this_surface.moveBy([scan_pos 0 0]);
                this_circle = circle.move([scan_pos 0 0]);
            case 'y'
                this_surface.moveBy([scan_pos*tand(init_angle) -scan_pos 0]);
                this_circle = circle.move([scan_pos*tand(init_angle) -scan_pos 0]);
            case 'z'
                this_surface.moveBy([0 0 scan_pos]);
                this_circle = circle.move([0 0 scan_pos]);
            otherwise
                error('Specify a correct direction for the line scan')
        end
        
        % Direct beam
        [~, killed, numScattersRay] = switch_plate('plate_represent', ...
            plate_represent, 'sample', this_surface, 'max_scatter', maxScatter, ...
            'pinhole_surface', pinhole_surface, 'thePlate', thePlate, ...
            'sphere', sphere, 'circle', this_circle, 'ray_model', ...
            ray_model, 'which_beam', direct_beam.source_model, 'beam', direct_beam);

        % Effuse beam
        [~, effuseKilled, numScattersEffuse] = switch_plate('plate_represent', ...
            plate_represent, 'sample', this_surface, 'max_scatter', maxScatter, ...
            'pinhole_surface', pinhole_surface, 'thePlate', thePlate, ...
            'sphere', sphere, 'circle', this_circle, 'ray_model', ...
            ray_model, 'which_beam', 'Effuse', 'beam', effuse_beam);

        % Update the progress bar if we are working in the MATLAB GUI.
        if progressBar
            if ~isOctave
                ppm.increment();
            elseif isOctave
                waitbar(i_/n_pixels, h);
            end
        end
        
        % Save the data for this iteration
        cntr_effuse_single(:,i_) = numScattersEffuse(1);
        counter_effuse_multiple(:,i_) = sum(numScattersEffuse(2:end));
        killed_effuse(i_) = effuseKilled;
        counters(:,:,i_) = numScattersRay;
        num_killed(i_) = killed;

        % Delete the surface object for this iteration
        if ~isOctave
            delete(this_surface);
        end
    end

    % Close the parallel pool
    % current_pool = gcp('nocreate');
    % delete(current_pool);

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
    if strcmp(scan_inputs.direction_1D, 'y')
        scan_range = scan_inputs.range1D + dist_to_sample;
    end
    line_scan_info = LineInfo(scan_inputs.direction_1D, scan_range, dist_to_sample, counters, num_killed, ...
        scan_inputs.raster_movment1D, direct_beam.n, t, t_estimate, cntr_effuse_single, ...
        counter_effuse_multiple, killed_effuse, direct_beam);

    if progressBar && ~isOctave
        delete(ppm);
    end
    
    if progressBar && make_plots
        line_scan_info.producePlots(thePath);
    end
end
