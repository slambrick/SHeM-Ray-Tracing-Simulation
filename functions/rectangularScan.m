% rectangularScan.m
%
% Copyright (c) 2018-19, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.
%
% Generates a 2d simulated image. Utilises a parfor loop to parellise the code.
%
% TODO: name value pair inputs
%
% Calling syntax: 
%
% INPUTS:
%
% OUTPUTS:
%  square_scan_info - An object of class RectangleInfo that contains all
%                     the results and information about the simulation
function square_scan_info = rectangularScan(sample_surface, xrange, zrange, ...
        direct_beam, raster_movement_x, raster_movement_z, maxScatter, ...
        pinhole_surface, effuse_beam, dist_to_sample, ...
        sphere, thePath, pinhole_model, thePlate, apertureAbstract, ray_model, ...
        n_detector)
    
    % The sample positions
    sample_xs = xrange(1):raster_movement_x:xrange(2);
    sample_zs = zrange(1):raster_movement_z:zrange(2);
    nx_pixels = length(sample_xs);
    nz_pixels = length(sample_zs);
    
    % Move the sample to the corner of the positions
    sample_surface.moveBy([xrange(1) 0 zrange(1)]);
    sphere.c(1) = sphere.c(1) + xrange(1);
    sphere.c(3) = sphere.c(3) + zrange(1);
    
    % Create the variables for output data
    counters = zeros(maxScatter, n_detector, nz_pixels, nx_pixels);
    effuse_counters = zeros(n_detector, nz_pixels, nx_pixels);
    num_killed = zeros(nz_pixels, nx_pixels);
    
    % Produce a time estimage for the simulation and print it out. This is
    % nessacerily a rough estimate.
    N_pixels = nx_pixels*nz_pixels;
    
    % Estimate of the time for the simulation
    t_estimate = time_estimate('n_rays', direct_beam{1}, 'n_effuse', ...
        effuse_beam{1}, 'sample_surface', sample_surface, 'n_pixels', N_pixels, ...
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
        ppm = ParforProgMon('Simulation progress: ', N_pixels);
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
    parfor i_=1:N_pixels
        
        % The x and z pixels we are on
        z_pix = mod(i_, nz_pixels) - 1;
        if mod(i_, nz_pixels) == 0
            x_pix = floor(i_/nz_pixels) - 1;
        else
            x_pix = floor(i_/nz_pixels);
        end
        
        if z_pix == -1
            z_pix = nz_pixels - 1;
        end
        
        % Place the sample into the right position for this pixel
        this_surface = copy(sample_surface);
        this_surface.moveBy([raster_movement_x*x_pix, 0, raster_movement_z*z_pix]);
        this_sphere = sphere;
        this_sphere.c(1) = this_sphere.c(1) + raster_movement_x*x_pix;
        this_sphere.c(3) = this_sphere.c(3) + raster_movement_z*z_pix;
        
        % Direct beam
        [~, killed, numScattersRay] = switch_plate('plate_represent', ...
            plate_represent, 'sample', this_surface, 'maxScatter', maxScatter, ...
            'pinhole_surface', pinhole_surface, 'thePlate', thePlate, ...
            'dist', dist_to_sample, 'sphere', this_sphere, 'ray_model', ...
            ray_model, 'which_beam', direct_beam{5}, 'beam', direct_beam);
        
        % Effuse beam
        [effuse_cntr, ~, ~] = switch_plate('plate_represent', ...
            plate_represent, 'sample', this_surface, 'maxScatter', maxScatter, ...
            'pinhole_surface', pinhole_surface, 'thePlate', thePlate, ...
            'dist', dist_to_sample, 'sphere', this_sphere, 'ray_model', ...
            ray_model, 'which_beam', 'Effuse', 'beam', effuse_beam);
        
        % Update the progress bar if we are working in the MATLAB GUI.
        if progressBar && ~isOctave
            ppm.increment();
        elseif isOctave
            waitbar(i_/N_pixels, h);
        end
        
        counters(:,:,i_) = numScattersRay;
        num_killed(i_) = killed;
        effuse_counters(:,i_) = effuse_cntr';
        
        % Delete the surface object for this iteration
        if ~isOctave
            delete(this_surface);
        end
    end
    
    % Close the parallel pool
    if ~isOctave
        current_pool = gcp('nocreate');
        delete(current_pool);
    end
    
    t = toc;
    
    % Actual time taken
    fprintf('Actual time taken: %f s\n', t);
    hr = floor(t/(60^2));
    mins = floor((t - hr*60*60)/60);
    if mins == 60
        hr = hr + 1;
        mins = 0;
    end
    fprintf('Which is: %i hr %2i mins\n\n', hr, mins);
    
    % Generate output square scan class
    square_scan_info = RectangleInfo(counters, num_killed, sample_surface, ...
        xrange, zrange, raster_movement_x, raster_movement_z, direct_beam{1}, ...
        effuse_beam{1}, t, t_estimate, effuse_counters, n_detector, maxScatter);
    
    % Draw and save images
    % All the images are saved giving maximum contrast in the images:
    %  black - pixel with fewest counts
    %  white - pixel with the most counts
    % Only draws them if there is a GUI.
    if progressBar
        %square_scan_info.produceImages(thePath);
    end
end

