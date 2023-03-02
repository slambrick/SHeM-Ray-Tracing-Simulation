% Copyright (c) 2018-22, Sam Lambrick, Aleksander Radic.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.

%close all
clear
%clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Paths to functions
loadpath

%% Read and define parameters
% Parameters are read and parsed in a seperate script called by this one
define_parameters

%% Check some of the inputs.
% This is by no means exhaustive. Most checks are for obvious errors.
% TODO: update/fix
% NOTE: this is not used and woefully out of date
if false
    if dist_to_sample < 0.1
        error(['Cannot guarentee good results when the sample is very close to' ...
            'the pinhole plate.']);
    end

    if strcmp(typeScan, 'line') && (range1D(1) > range1D(2))
        error('Impossible range for a line scan specified.');
    end

    if strcmp(typeScan, 'line') && strcmp(Direction, 'y')
        if (dist_to_sample + range1D(1)) < 0.1
            error('y line scan goes too close or behind the pinhole plate.');
        end
    end

    if strcmp(typeScan, 'rectangular') && ...
            ((xrange(1) > xrange(2)) || (zrange(1) > zrange(2)))
        error('Impossible range for a 2D scan specified.')
    end

    if (~strcmp(sample_type, 'flat') && ~strcmp(sample_type, 'sphere') ...
            && (~strcmp(sample_type, 'custom')))
        error('Specify a correct type of sample.')
    end

    if ((strcmp(sample_type, 'sphere') || strcmp(sample_type, 'flate')) && ...
            square_size <=0)
        error('The size of the flat square must be positive and non-zero.');
    end

    if scale <= 0
        error('The scaling of the sample model must be positive and non-zero.')
    end

    if (~strcmp(pinhole_model, 'stl') && ~strcmp(pinhole_model, 'circle') && ...
            ~strcmp(pinhole_model, 'aperture') && ~strcmp(pinhole_model, 'new') ...
            && ~strcmp(pinhole_model, 'abstract'))
        error('Specify a correct model of pinhole .');
    end
end

%% Path for simulation results

% Tha path to save the simulation results to
thePath = simulationDir(directory_label);

if ~exist(thePath, 'dir')
    mkdir(thePath)
end
copyfile(param_fname, thePath)


% Are we running in GNU Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Add the required libraries if we are running in octave
% TODO: check the compatibility
if isOctave
    pkg load statistics;
    pkg load image;
end

addpath(thePath);

%% Sample import and plotting

% A struct to represent the sphere and circle parts of the sample
sphere = Sphere(1, material, sphere_cs, sphere_rs);
circle = Circle(1, material, circle_c, circle_r, circle_n);

% Importing the sample as a TriagSurface object.
[sample_surface, sphere, circle, sample_description] = sample_import(sample_inputs, sphere, ...
    circle, pinhole_plate_inputs.working_dist, dontMeddle, square_size, init_angle);
sample_inputs.sample_descrition = sample_description;

if strcmp(typeScan, 'line')
    sample_surface.moveBy(init_displacement);
end

% Do any extra manipulation of the sample here

% Translation of the sphere
if false
    % This is for indenting the sphere into the surface for the 3D
    % reconstruction B-SHeM test.
    sphere.centre = sphere.centre + [0;-1.36;0];
end

if false
    % Tilt required to explain the problem with displacement of the diffraction
    % p[attern
    sample_surface.rotateGeneral('x', -2.1);
    sample_surface.rotateGeneral('z', -3.8);
end

if false
    sample_surface.rotateGeneral('y', -60+90+20);
    sample_surface.rotateGeneral('y', 180);
    sample_surface.moveBy([0.05, 0, 0])
    %sample_surface.moveBy([0, 0.525, 0]);
end

% Specifically for the simulation of the LiF diffrtaction with multiscat
% peak intensities
if false
    % Flip the z coordinates of the lattice vectors because of a difference
    % in the axes used in multiscat
    sample_surface.lattice(:,3) = sample_surface.lattice(:,3);
    sample_surface.lattice(:,6) = -sample_surface.lattice(:,6);
end
%sample_surface.rotateGeneral('y', 45);

% Plot the sample surface in 3D space, if we are using a graphical window
% TODO: put in a seperate
if feature('ShowFigureWindows')
    if ~strcmp(typeScan, 'single_pixel') && ~strcmp(sample_inputs.sample_type, 'circle')
        sample_surface.patchPlot(true);
        ylim([-dist_to_sample - 0.2, -dist_to_sample + 0.2]);
    end

    if strcmp(typeScan, 'rectangular') || strcmp(typeScan, 'rotations') ||...
            strcmp(typeScan, 'multiple_rectangular')
        xlim(xrange);
        zlim(zrange);
        ylim(zrange - dist_to_sample);
    elseif strcmp(typeScan, 'line') && ~strcmp(Direction, 'y')
        xlim(range1D);
        zlim(range1D);
        ylim(range1D - dist_to_sample);
    elseif strcmp(typeScan, 'line') && strcmp(Direction, 'y')
        xlim([-0.2, 0.2]);
        zlim([-0.2, 0.2]);
        ylim([-0.2, 0.2] - dist_to_sample);

    end

    if ~strcmp(typeScan, 'single pixel')
        print(fullfile(thePath, 'sample_closeUp.eps'), '-depsc');
    end
end

%sample_surface.normalise_lattice();

%% Pinhole plate import and plotting
[pinhole_surface, thePlate, pinhole_model] = pinhole_import(...
    pinhole_plate_inputs, sample_surface, defMaterial);

if thePlate.in_use
    thePlate.plot_detectors;
end
%% Compile the mex files

mexCompile(recompile);


%% Performing the simulation
% TODO: update
disp(sample_surface)
switch typeScan
    case 'rectangular'
        % For a rectangular scan
        if init_angle_pattern
            raster_pattern = generate_raster_pattern('raster_movment2D', ...
                [raster_movment2D_x, raster_movment2D_z], 'xrange', xrange, ...
                'zrange', zrange, 'init_angle', init_angle);
        else
            raster_pattern = generate_raster_pattern('raster_movment2D', ...
                [raster_movment2D_x, raster_movment2D_z], 'xrange', xrange, ...
                'zrange', zrange);
        end
        simulationData = rectangularScan('sample_surface', sample_surface, ...
            'raster_pattern', raster_pattern,'direct_beam', direct_beam, ...
            'max_scatter', max_scatter,      'pinhole_surface', pinhole_surface, ...
            'effuse_beam', effuse_beam,      'dist_to_sample', dist_to_sample, ...
            'circle', circle, ...
            'sphere', sphere,                'thePath', thePath, ...
            'pinhole_model', pinhole_model,  'thePlate', thePlate, ...
            'ray_model', ray_model,          'n_detector', n_detectors);
    case 'multiple_rectangular'
        % TODO: check this works and then make it work with the new parameter
        % specification file
        simulationData = {};
        
        % find y positions
        ys = range_y(1):raster_movement_y:range_y(2);
        ny = length(ys);

        % progress bar
        h = waitbar(0, 'Proportion of simulations performed', 'Name', 'Ray tracing progress');

        for i_ = 1:ny
            y_displacement = ys(i_);

            % Move the sample in y by given amount
            surface_copy = copy(sample_surface);
            surface_copy.moveBy([y_displacement, -y_displacement, 0]);
            sphere.centre = sphere.centre + [y_displacement, -y_displacement, 0];
            y_distance = dist_to_sample + y_displacement;
            
            % Create the raster pattern at this z
            if init_angle_pattern
                raster_pattern = generate_raster_pattern('raster_movment2D', ...
                    [raster_movment2D_x, raster_movment2D_z], 'xrange', xrange, ...
                    'zrange', zrange, 'init_angle', init_angle);
            else
                raster_pattern = generate_raster_pattern('raster_movment2D', ...
                    [raster_movment2D_x, raster_movment2D_z], 'xrange', xrange, ...
                    'zrange', zrange);
            end
            
            % Create a directory for this scan
            subPath = sprintf('z%.4fmm', y_displacement + dist_to_sample);
            subPath = fullfile(thePath, subPath);
            if ~exist(subPath, 'dir')
                mkdir(subPath)
            end

            % Run the simulation at that z
            simulationData{i_} = rectangularScan('sample_surface', surface_copy, ...
                'raster_pattern', raster_pattern,'direct_beam', direct_beam, ...
                'max_scatter', max_scatter,      'pinhole_surface', pinhole_surface, ...
                'effuse_beam', effuse_beam,      'dist_to_sample', y_distance, ...
                'sphere', sphere,                'thePath', subPath, ...
                'circle', circle, ...
                'pinhole_model', pinhole_model,  'thePlate', thePlate, ...
                'ray_model', ray_model,          'n_detector', n_detectors); %#ok<SAGROW>
            waitbar(i_/ny, h);
        end

        close(h);
        delete(h);
    case 'line'
        % For a line scan
        % TODO: update with the new lower level functions
        simulationData = lineScan('sample_surface', sample_surface, ...
            'scan_inputs', scan_inputs,     'direct_beam', direct_beam, ...
            'max_scatter', max_scatter,     'pinhole_surface', pinhole_surface, ...
            'effuse_beam', effuse_beam,     'dist_to_sample', dist_to_sample, ...
            'sphere', sphere,               'thePath', thePath, ...
            'circle', circle, ...
            'pinhole_model', pinhole_model, 'thePlate', thePlate, ...
            'ray_model', ray_model,         'n_detector', n_detectors, ...
            'make_plots', true,             'init_angle', init_angle);
    case 'single pixel'
        % For a single pixel
        % TODO: update with the new lower level functions
        simulationData = singlePixel(sample_surface, direct_beam, ...
            max_scatter, pinhole_surface, effuse_beam, dist_to_sample, sphere, ...
            thePath, save_to_text, pinhole_model, thePlate, aperture_abstract);
    case 'rotations'
        % Perform multiple scans while rotating the sample in between.
        % TODO: add the ability to rotate about a different axis
        simulationData = {};
        h = waitbar(0, 'Proportion of simulations performed', 'Name', 'Ray tracing progress');
        N = length(rot_angles);
        sphere_centre = sphere.centre;
        
        % Loop through the rotations
        for i_=1:N
            s_surface = copy(sample_surface);
            s_surface.rotateGeneral('y', rot_angles(i_));
            
            % Rotate the centre of the sphere
            theta = rot_angles(i_)*pi/180;
            s = sin(theta);
            c = cos(theta);
            R = [c, 0, s; 0, 1, 0; -s, 0, c];
            sphere.centre = (R*sphere_centre)';
            
            subPath = [thePath '/rotation' num2str(rot_angles(i_))];
            if ~exist(subPath, 'dir')
                mkdir(subPath)
            end
            
            tmp_angle = init_angle;
            if init_angle_pattern
                init_angle = 0;
            end
            
            switch scan_pattern
                case 'rotations'
                    raster_pattern = generate_raster_pattern('raster_movment2D', ...
                        [raster_movment2D_x, raster_movment2D_z], 'xrange', xrange, ...
                        'zrange', zrange, 'rot_angle', rot_angles(i_), 'init_angle', init_angle);
                case 'regular'
                    raster_pattern = generate_raster_pattern('raster_movment2D', ...
                        [raster_movment2D_x, raster_movment2D_z], 'xrange', xrange, ...
                        'zrange', zrange, 'init_angle', init_angle);
                otherwise
                    error('Please specify an existing scan pattern')
            end
            init_angle = tmp_angle;
            
            simulationData{i_} = rectangularScan('sample_surface', s_surface, ...
                'raster_pattern', raster_pattern,'direct_beam', direct_beam, ...
                'max_scatter', max_scatter,      'pinhole_surface', pinhole_surface, ...
                'effuse_beam', effuse_beam,      'dist_to_sample', dist_to_sample, ...
                'sphere', sphere,                'thePath', subPath, ...
                'circle', circle, ...
                'pinhole_model', pinhole_model,  'thePlate', thePlate, ...
                'ray_model', ray_model,          'n_detector', n_detectors); %#ok<SAGROW>
            
            waitbar(i_/N, h);
            
            % Save all data to a .mat file
            if output_data
                save([thePath '/' data_fname], 'simulationData', 'sample_inputs', ...
                    'direct_beam', 'effuse_beam', 'pinhole_plate_inputs', 'scan_inputs');
            end

        end
        
        if strcmp(scan_pattern, 'rotations')
            re_rotate_images(simulationData, thePath, rot_angles)
        end
        close(h);
        delete(h);
    case 'line_rotations'
        % Perform multiple scans while rotating the sample in between.
        simulationData = {};
        h = waitbar(0, 'Proportion of simulations performed');
        N = length(rot_angles);
        sphere_centre = sphere.centre;

        % Loop through the rotations
        for i_=1:N
            %close all % <- there are more elegant ways of doing this...
            s_surface = copy(sample_surface);
            s_surface.rotateGeneral('y', -rot_angles(i_));
            
            % TODO: change the material of the sample to have the correct
            % parameters for the scattering distribution
            if true
                th = 0;
                if rot_angles(i_) > 90 && rot_angles(i_) <= 180
                    s_surface.rotateLattice('y', 90);
                    th = rot_angles(i_) - 90;
                elseif rot_angles(i_) > 180 && rot_angles(i_) <= 270
                    s_surface.rotateLattice('y', 180);
                    th = rot_angles(i_) - 180;
                elseif rot_angles(i_) > 270 && rot_angles(i_) <= 360
                    s_surface.rotateLattice('y', 270);
                    th = rot_angles(i_) - 270;
                end
                if rot_angles(i_) > 90
                    for j_=1:length(s_surface.nTriag)
                        s_surface.compositions{j_} = ['diffractive_' num2str(th) 'deg'];
                    end
                end
            end

            % Rotate the centre of the sphere
            theta = rot_angles(i_)*pi/180;
            s = sin(theta);
            c = cos(theta);
            R = [c, 0, s; 0, 1, 0; -s, 0, c];
            sphere.centre = (R*sphere_centre)';
            %s_surface.patchPlot(true);

            subPath = [thePath '/rotation' num2str(rot_angles(i_))];
            if ~exist(subPath, 'dir')
                mkdir(subPath)
            end
            %disp(s_surface.lattice);
            simulationData{i_} = lineScan('sample_surface', s_surface, ...
                'scan_inputs', scan_inputs,     'direct_beam', direct_beam, ...
                'max_scatter', max_scatter,     'pinhole_surface', pinhole_surface, ...
                'effuse_beam', effuse_beam,     'dist_to_sample', dist_to_sample, ...
                'sphere', sphere,               'thePath', subPath, ...
                'circle', circle, ...
                'pinhole_model', pinhole_model, 'thePlate', thePlate, ...
                'ray_model', ray_model,         'n_detector', n_detectors, ...
                'make_plots', false,            'init_angle', init_angle); %#ok<SAGROW>
            
            waitbar(i_/N, h);
        end
        close(h);
        delete(h);
    otherwise
        error(['Need to specify a valid type of scan: "line", ', ...
               '"rectangular", "single pixel"']);
end

%% Output data about simulation to files

% Save all data to a .mat file
if output_data
    save(fullfile(thePath, data_fname), 'simulationData', 'sample_inputs', ...
        'direct_beam', 'effuse_beam', 'pinhole_plate_inputs', 'scan_inputs');
end

% Save formatted data to a .mat file, does not include all the parameters
% but does include the core outputs.
if contains(typeScan, {'rotations', 'rectangular', 'line', 'line_rotations'})
    if output_data && (strcmp(typeScan, 'rotations') || strcmp(typeScan, 'line_rotations'))
        formatOutputRotation(simulationData, thePath, rot_angles);
    elseif output_data
        simulationData.formatOutput(thePath);
    end
end

% Save main data to a text file
textFname = 'data_for_plotting.csv';
if save_to_text && strcmp(typeScan, 'rotations')
    for i_=1:length(rot_angles)
        % TODO: bug here
        currentFname = [textFname(1:end-4) num2str(rot_angles(i_)) '.csv'];
        simulationData{i_}.saveText([thePath '/' currentFname num2str(rot_angles(i_))]);
    end
elseif save_to_text && ~(strcmp(typeScan, 'line_rotations') || ...
        strcmp(typeScan, 'rotations') || strcmp(typeScan, 'multiple_rectangular'))
    simulationData.saveText([thePath '/' textFname]);
end

% Save the parameters to a text file
%if saveParams
%    saveParam(sample_inputs, direct_beam, effuse_beam, ...
%        pinhole_plate_inputs, scan_inputs);
%end

