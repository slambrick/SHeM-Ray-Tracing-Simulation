% Copyright (c) 2018-21, Sam Lambrick, Aleksander Radic.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.

%close all
clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Paths to functions
loadpath

%% Read parameters from text file
% TODO: move parameter reading into a seperate file
param_fname = 'ray_tracing_parameters.txt';
param_list = read_parameters(param_fname);

% Set up virtual microscope
working_dist = str2double(param_list{1});
init_angle = str2double(param_list{2});
typeScan = strtrim(param_list{3});
n_detectors = str2double(param_list{4});
aperture_axes = parse_list_input(param_list{5});
aperture_c = parse_list_input(param_list{6});
rot_angles = parse_rotations(param_list{7});
pinhole_model = parse_pinhole(param_list{8});

% Set up source
n_rays = str2double(param_list{9});
pinhole_r = str2double(param_list{10});
source_model = strtrim(param_list{11});
theta_max = str2double(param_list{12});
sigma_source = str2double(param_list{13});
if ~parse_yes_no(param_list{14})
    effuse_size = 0;
else
    effuse_size = str2double(param_list{15});
end

% Set up sample
sample_type = strtrim(param_list{16});
material = parse_scattering(strtrim(param_list{17}), str2double(param_list{18}), ...
    str2double(param_list{19}));
sample_description = param_list{20};
dist_to_sample = str2double(param_list{21});
sphere_r = str2double(param_list{22});
square_size = str2double(param_list{23});
sample_fname = strtrim(param_list{24});
dontMeddle = parse_yes_no(param_list{25});

% Set up scan
pixel_seperation = str2double(param_list{26});
range_x = str2double(param_list{27});
range_z = str2double(param_list{28});
init_angle_pattern = ~parse_yes_no(param_list{29});

% 1D scan parameters
scan_direction = strtrim(param_list{30});
range_1D = [str2double(param_list{31}), str2double(param_list{32})];
res_1D = str2double(param_list{33});

% Other parameters
directory_label = strtrim(param_list{34});
recompile = parse_yes_no(param_list{35});

%% Generate parameters from the inputs 

pinhole_c = [-working_dist*tand(init_angle), 0, 0];
n_effuse = n_rays*effuse_size;
raster_movment2D_x = pixel_seperation;
raster_movment2D_z = pixel_seperation;
xrange = [-range_x/2, range_x/2];
zrange = [-range_z/2, range_z/2];
sphere_c = [0, -dist_to_sample + sphere_r, 0];

%% Define remaining parameters

% The maximum number of sample scatters per ray. There is a hard-coded total
% maximum number of scattering events of 1000 (sample and pinhole plate). Making
% this uneccaserily large will increase the memory requirments of the
% simulation.
max_scatter = 20;

% If rotations are present the scan pattern can be regular or be adjusted to
% match the rotation of the sample
scan_pattern = 'regular';

% Do we want to generate rays in Matlab (more flexibility, more output options)
% or in C (much lower memory requirments and slightly faster), 'C' or 'MATLAB'
% In general stick to 'C' unless your own source model is being used
ray_model = 'C';

% Exponant of the cosine in the effuse beam model
cosine_n = 1;

% In the case of the predefined CAD model, specify the accuraccy of the
% triangulation, 'low', 'medium', or 'high' (use 'low').
plate_accuracy = 'low';

% In the case of 'circle', specify the radius of the circle (mm).
circle_plate_r = 4;

% Should a flat pinhole plate be modelled (with 'N circle'). not including may
% speed up the simulation but won't model the effuse and multiple scattering
% backgrounds properly.
plate_represent = 0;

% In the case of 'abstract', specify the two angles of the location of the
% detector aperture and the half cone angle of its extent. Note that the
% aperture can only be placed in the hemisphere facing the sample. All
% angles in degrees.
% TODO: this
aperture_theta = 0;
aperture_phi = 0;
aperture_half_cone = 15;


%% Parameters for multiple rectangular scans
raster_movement_y = 1.4;   % increment between 2 scans
range_y = [0    1.4];      % range of y positions relative to 'dist_to_sample'
% For line scans in the y-direction be careful that the sample doesn't go
% behind the pinhole plate.
init_displacement = [0, 0, 0];  % initial position of sample from 'centred'
raster_movment1D = 0.02;        % movement increment
range1D = [-1 4];               % range
Direction = 'y';                % 'x', 'y' or 'z' - along which direction to move

% make the model 10 times smaller (Inventor exports in cm by default...).
scale = 2;
%  'strips' - Two parallel series of strips with varying parameters
% make the model 10 times larger (Inventor exports in cm by default...).
scale = 1;

% A string giving a brief description of the sample, for use with
% sample_type = 'custom'
sample_description = 'Sample with series of diffractive peaks.';

%% Parameters of the default material, to use for faces where no material is specified
defMaterial.function = 'cosine';
defMaterial.params = 0;
defMaterial.color = [0.8 0.8 1.0];


%% Output and plotting parameters

% Which figures to plot
% The starting positions of the rays and the number of rays at each point
plot_density = false;

% Save simulation data to file
output_data = true;
data_fname = 'scatteringData.mat';
saveParams = true;
paramsFile = 'scatteringParameters.txt';
% Applies only to line scans, saves in results_path with name 'data_for_plotting.csv'
save_to_text = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create input structs to hold input data

% TODO: convert from structs into classes??
direct_beam.n = n_rays;
direct_beam.pinhole_c = pinhole_c;
direct_beam.pinhole_r = pinhole_r;
direct_beam.theta_max = theta_max;
direct_beam.source_model = source_model;
direct_beam.init_angle = init_angle;
direct_beam.sigma_source = sigma_source;

effuse_beam.n = n_effuse;
effuse_beam.pinhole_c = pinhole_c;
effuse_beam.pinhole_r = pinhole_r;
effuse_beam.cosine_n = cosine_n;

% NOTE: I think this has to remain in the main script, however it could be
% placed earlier and then only a few parameter structs would need to be passed
% around... ??
scan_inputs.type_scan = typeScan;
scan_inputs.max_scatter = max_scatter;
switch typeScan
    case 'multiple_rectangular'
        scan_inputs.rotationAngles = 0;
        scan_inputs.raster_movment2D_x = raster_movment2D_x;
        scan_inputs.raster_movment2D_z = raster_movment2D_z;
        scan_inputs.xrange = xrange;
        scan_inputs.zrange = zrange;
        scan_inputs.raster_movment1D = raster_movement_y;
        scan_inputs.range1D = range_y;
        scan_inputs.direction_1D = 'y';
    case 'rotations'
        scan_inputs.rotationAngles = rot_angles;
        scan_inputs.raster_movment2D_x = raster_movment2D_x;
        scan_inputs.raster_movment2D_z = raster_movment2D_z;
        scan_inputs.xrange = xrange;
        scan_inputs.zrange = zrange;
        scan_inputs.raster_movment1D = NaN;
        scan_inputs.range1D = NaN;
        scan_inputs.direction_1D = NaN;
    case 'rectangular'
        scan_inputs.rotationAngles = 0;
        scan_inputs.raster_movment2D_x = raster_movment2D_x;
        scan_inputs.raster_movment2D_z = raster_movment2D_z;
        scan_inputs.xrange = xrange;
        scan_inputs.zrange = zrange;
        scan_inputs.raster_movment1D = NaN;
        scan_inputs.range1D = NaN;
        scan_inputs.direction_1D = NaN;
    case 'line'
        scan_inputs.rotationAngles = 0;
        scan_inputs.raster_movment2D_x = NaN;
        scan_inputs.raster_movment2D_z = NaN;
        scan_inputs.xrange = NaN;
        scan_inputs.zrange = NaN;
        scan_inputs.raster_movment1D = res_1D;
        scan_inputs.range1D = range_1D;
        scan_inputs.direction_1D = scan_direction;
    case 'line_rotations'
        scan_inputs.rotationAngles = rot_angles;
        scan_inputs.raster_movment2D_x = NaN;
        scan_inputs.raster_movment2D_z = NaN;
        scan_inputs.xrange = NaN;
        scan_inputs.zrange = NaN;
        scan_inputs.raster_movment1D = res_1D;
        scan_inputs.range1D = range_1D;
        scan_inputs.direction_1D = scan_direction;
    case 'single pixel'
        scan_inputs.rotationAngles = 0;
        scan_inputs.raster_movment2D_x = NaN;
        scan_inputs.raster_movment2D_z = NaN;
        scan_inputs.xrange = NaN;
        scan_inputs.zrange = NaN;
        scan_inputs.raster_movment1D = NaN;
        scan_inputs.range1D = NaN;
        scan_inputs.direction_1D = NaN;
end

sample_inputs.sample_type = sample_type;
sample_inputs.material = material; % Only use this for stl or homebuilt surfaces
sample_inputs.dist_to_sample = dist_to_sample;
sample_inputs.sphere = sphere;
sample_inputs.sample_description = sample_description;
if strcmp(sample_type, 'rectangular')
    sample_inputs.square_size = square_size;
else
    sample_inputs.square_size = NaN;
end
if strcmp(sample_type, 'custom')
    sample_inputs.sample_fname = sample_fname;
else
    sample_inputs.sample_fname = NaN;
end
sample_inputs.scale = scale;

pinhole_plate_inputs.pinhole_model = pinhole_model;
pinhole_plate_inputs.working_dist = working_dist;
switch pinhole_model
    case {'stl', 'new', 'angular'}
        pinhole_plate_inputs.plate_accuracy = plate_accuracy;
        pinhole_plate_inputs.n_detectors = 1;
        pinhole_plate_inputs.plate_represent = 1;
        pinhole_plate_inputs.aperture_axes = NaN;
        pinhole_plate_inputs.aperture_c = NaN;
        pinhole_plate_inputs.circle_plate_r = NaN;
        pinhole_plate_inputs.aperture_theta = NaN;
        pinhole_plate_inputs.aperture_phi = NaN;
        pinhole_plate_inputs.aperture_half_cone = NaN;
    case 'circle'
        % NOTE: this is depreciated, probably remove
        pinhole_plate_inputs.plate_accuracy = NaN;
        pinhole_plate_inputs.n_detectors = 1;
        pinhole_plate_inputs.plate_represent = plate_represent;
        pinhole_plate_inputs.aperture_axes = aperture_axes;
        pinhole_plate_inputs.aperture_c = aperture_c;
        pinhole_plate_inputs.circle_plate_r = circle_plate_r;
        pinhole_plate_inputs.aperture_theta = NaN;
        pinhole_plate_inputs.aperture_phi = NaN;
        pinhole_plate_inputs.aperture_half_cone = NaN;
    case 'N circle'
        pinhole_plate_inputs.plate_accuracy = NaN;
        pinhole_plate_inputs.n_detectors = n_detectors;
        pinhole_plate_inputs.plate_represent = plate_represent;
        pinhole_plate_inputs.aperture_axes = aperture_axes;
        pinhole_plate_inputs.aperture_c = aperture_c;
        pinhole_plate_inputs.circle_plate_r = circle_plate_r;
        pinhole_plate_inputs.aperture_theta = NaN;
        pinhole_plate_inputs.aperture_phi = NaN;
        pinhole_plate_inputs.aperture_half_cone = NaN;
    case 'abstract'
        % TODO: make the 'abstract' simulations work
        pinhole_plate_inputs.n_detectors = n_detectors;
        pinhole_plate_inputs.plate_represent = 0;
        pinhole_plate_inputs.circle_plate_r = NaN;
        pinhole_plate_inputs.aperture_theta = NaN;
        pinhole_plate_inputs.aperture_phi = NaN;
        pinhole_plate_inputs.aperture_half_cone = NaN;
end

% Input structs are:
%  - pinhole_plate_inputs
%  - sample_inputs
%  - scan_inputs
%  - direct_beam
%  - effuse_beam
% I think they should be converted into classes
% There can then be an overall struct/object of the inputs
inputs.sample = sample_inputs;
inputs.pinhole = pinhole_plate_inputs;
inputs.scan = scan_inputs;
inputs.effuse_beam = effuse_beam;
inputs.direct_beam = direct_beam;
% TODO: use inputs instead of the seperate structs

%% Check some of the inputs.
% This is by no means exhaustive. Most checks are for obvious errors.
% TODO: update/fix
if false
    if dist_to_sample < 0.1
        error(['Cannot guarentee good results when the sample is very close to' ...
            'the pinhole plate.']);
    end

    if sphere_r*2 > dist_to_sample
        error('The sphere is too big for the pinhole plate-sample distance.');
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

    if (strcmp(sample_type, 'sphere') && sphere_r <= 0)
        error('The sphere must have a positive non-zero radius.')
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

% A struct to represent the sphere
sphere = Sphere(1, defMaterial, sphere_c, sphere_r);

% Importing the sample as a TriagSurface object.
[sample_surface, sphere, sample_description] = sample_import(sample_inputs, sphere, ...
    pinhole_plate_inputs.working_dist, dontMeddle, square_size);
sample_inputs.sample_descrition = sample_description;

if strcmp(typeScan, 'line')
    sample_surface.moveBy(init_displacement);
end

% Do any extra manipulation of the sample here
if true
    sample_surface.rotateGeneral('x', 4);
    sample_surface.rotateGeneral('z', 4);
end

% Plot the sample surface in 3D space, if we are using a graphical window
% TODO: put in a seperate
if feature('ShowFigureWindows')
    if ~strcmp(typeScan, 'single_pixel')
        sample_surface.patchPlot(true);
        ylim([-dist_to_sample - 0.2, -dist_to_sample + 0.2]);
    end

    if strcmp(typeScan, 'rectangular') || strcmp(typeScan, 'rotations') ||...
            strcmp(typeScan, 'multiple_rectangular')
        xlim([-xrange(2) -xrange(1)]);
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
[pinhole_surface, thePlate, aperture_abstract, pinhole_model] = pinhole_import(pinhole_plate_inputs, sample_surface);
%thePlate.backwall_represent = 1;
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
            'pinhole_model', pinhole_model, 'thePlate', thePlate, ...
            'ray_model', ray_model,         'n_detector', n_detectors, ...
            'make_plots', true);
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
            sphere.centre = (R*sphere_centre')';
            
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
            s_surface.rotateGeneral('y', rot_angles(i_));

            % Rotate the centre of the sphere
            theta = rot_angles(i_)*pi/180;
            s = sin(theta);
            c = cos(theta);
            R = [c, 0, s; 0, 1, 0; -s, 0, c];
            sphere.centre = (R*sphere_centre')';
            %s_surface.patchPlot(true);

            subPath = [thePath '/rotation' num2str(rot_angles(i_))];
            if ~exist(subPath, 'dir')
                mkdir(subPath)
            end
            disp(s_surface.lattice);
            simulationData{i_} = lineScan('sample_surface', s_surface, ...
                'scan_inputs', scan_inputs,     'direct_beam', direct_beam, ...
                'max_scatter', max_scatter,     'pinhole_surface', pinhole_surface, ...
                'effuse_beam', effuse_beam,     'dist_to_sample', dist_to_sample, ...
                'sphere', sphere,               'thePath', subPath, ...
                'pinhole_model', pinhole_model, 'thePlate', thePlate, ...
                'ray_model', ray_model,         'n_detector', n_detectors, ...
                'make_plots', false); %#ok<SAGROW>
            
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
        simulationData.saveText([thePath '/' currentFname]);
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

