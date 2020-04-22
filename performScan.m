% Copyright (c) 2018, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the
% GNU/GPL-3.0-or-later.

clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% General Parameters

% The maximum number of sample scatters per ray. There is a hard-coded total
% maximum number of scattering events of 1000 (sample and pinhole plate). Making
% this uneccaserily large will increase the memory requirments of the
% simulation.
maxScatter = 20;

% Type of scan 'line', 'rectangular', 'multiple_rectangular', 'rotations', or 'single pixel'
typeScan = 'multiple_rectangular';

% Recompile mex files?
% Required if using on a new computer or if changes to .c files have been made.
recompile = true;

%% Beam/source parameters %%
n_rays = 1e4;

% The inicidence angle in degrees
init_angle = 45;

% Geometry of pinhole
pinhole_c = 2.1*[-tand(init_angle), 0, 0];
pinhole_r = 0.6e-3;

% Number of rays to use and the width of the source

% skimmer radius over source - pinhole distance
theta_max = atan(0.01/100);

% Standard deviationo of the Gaussian model of the source
sigma_source = 0.1;

% Model for the virtual source to use, 'Gaussian' or 'Uniform'
source_model = 'Uniform';

% Put the information in a cell array to pass through functions
% TODO: use a struct rather than a cell array.
direct_beam.n = n_rays;
direct_beam.pinhole_c = pinhole_c;
direct_beam.pinhole_r = pinhole_r;
direct_beam.theta_max = theta_max;
direct_beam.source_model = source_model;
direct_beam.init_angle = init_angle;
direct_beam.sigma_source = sigma_source;

% Do we want to generate rays in Matlab (more flexibility, more output options)
% or in C (much lower memory requirments and slightly faster), 'C' or 'MATLAB'
ray_model = 'C';

%% Effusive beam parameters

% If a cosine is being used then specify the exponent of the cosine
cosine_n = 1;

% How large is the effusive beam (proportion of the size of the main beam). Set
% to zero if the effusive beam is not to be moddeled
effuse_size = 0;

% Information on the effuse beam
n_effuse = n_rays*effuse_size;
% TODO: use a struct rather than a cell array.
effuse_beam.n = n_effuse;
effuse_beam.pinhole_c = pinhole_c;
effuse_beam.pinhole_r = pinhole_r;
effuse_beam.cosine_n = cosine_n;

%% Pinhole plate parameters
% Specify how to model the pinhole plate:
%  'stl'       - Use the predefined CAD model of the pinhole plate (plate as it
%                is Feb 2018)
%  'circle'    - Use the detector aperture and model the pinhole plate as a
%                circle
%  'N circle'  - There are N circles in a plane as detector apertures
%  'new'       - CAD model of the pinhole plate in the new sample chamber (Apr
%                2019)
%  'new_micro' - TODO
%  'abstract'  - TODO
pinhole_model = 'N circle';

% In the case of the predefined CAD model, specify the accuraccy of the
% triangulation, 'low', 'medium', or 'high' (use 'low').
plate_accuracy = 'low';

% In the case of 'circle', specify the radius of the circle (mm).
circle_plate_r = 4;

% In the case of 'aperture' or 'circle' specify the axes of the aperture. Axis 1
% is along the beam direction ('x') and axis 2 is perpendicular to the beam
% direction ('z'). The aperture is always centred on the x-axis and is displaced
% by the specified amount.
n_detectors = 1;
aperture_axes = [1*sqrt(2)    1];
aperture_c = [2.1    0];
plate_represent = 0;

% In the case of 'abstract', specify the two angles of the location of the
% detector aperture and the half cone angle of its extent. Note that the
% aperture can only be placed in the hemisphere facing the sample. All
% angles in degrees.
aperture_theta = 0;
aperture_phi = 0;
aperture_half_cone = 15;

%% Parameters for a 2d scan
% Ususally the ranges should go from -x to x. Note that these limits are in the
% coordiante system of the final image - the x axis of the final image is the
% inverse of the simulation x axis.
raster_movment2D_x = 2e-3;
raster_movment2D_z = 2e-3;
xrange = [-0.15      0.15];
zrange = [-0.10      0.10];

%% Parameters for multiple rectangular scans
raster_movement_y = 0.1;   % increment between 2 scans
range_y = [-1   4];        % range of y positions relative to 'dist_to_sample'

%% Rotating parameters
% Parameters for multiple images while rotating the sample.
rot_angles = [0, 72, 144, 216, 288];

%% Parameters for a 1d scan
% For line scans in the y-direction be careful that the sample doesn't go
% behind the pinhole plate.
init_displacement = [0, 0, 0];  % initial position of sample from 'centred'
raster_movment1D = 0.02;        % movement increment
range1D = [-1 4];               % range
Direction = 'y';                % 'x', 'y' or 'z' - along which direction to move

%% Sample parameters

% What type of sample to use :
%  'flat'   - A flat square (need to specify square_size)
%  'strips' - Two parallel series of strips with varying parameters
%  'sphere' - An anlaytic sphere on a flat square surface (need to specify
%             square_size and sphere_r)
%  'photoStereo' - A test sample for photo stereo that includes part of a
%                  sphere along with input from a CAD file
%  'custom' - Uses a CAD model from file
%  'airy'   - TODO
%  'corrugation' - TODO
sample_type = 'custom';

% The sample file, include the full path
sample_fname = 'samples/lif_real.obj';

% Sample scaling, for if the CAD model had to be made at a larger scale. 10 will
% make the model 10 times larger (Inventor exports in cm by default...).
scale = 1;

% A string giving a brief description of the sample, for use with
% sample_type = 'custom'
sample_description = 'Realistic representation of the LiF sample.';

%% Parameters of the default material, to use for faces where no material is specified
defMaterial.function = 'cosine';
defMaterial.params = 0;
defMaterial.color = [0.8 0.8 1.0];

% The nominal working distance of the geometry
working_dist = 2.1;

% How close should the nearest point of the sample be to the pinhole plate, the
% defualt is 2.121 to maintain the 45o geometry. If an analytic sphere is being
% used then this is the distance between the flat surface the sphere sits on and
% the pinhole plate.
dist_to_sample = 2.1;

% The radius of the anayltic sphere (mm) (if it being included)
sphere_r = 0.1;

% Centre of the anayltic sphere (mm)
sphere_c = [0.05, -dist_to_sample - sphere_r*2/3, 0.05];

% If a flat sample is being used or if a sphere on a flat surface what is the
% length of the sides of the square.
square_size = 0.8;

%% Output and plotting parameters

% Where to save figures/data files
% All figures and output data will be saved to this directory.
directory_label = 'lif_real_mirrored';

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
        error('Specify a correct model of pinhole plate.');
    end
end

%% Paths to functions
addpath('import_3d/stlread', 'import_3d/objread', 'functions', ...
        'functions/interface_functions', 'classes', ...
        'mexFiles', 'DylanMuir-ParforProgMon-a80c9e9', 'functions/standard_samples', ...
        'surf2stl');

%% Path for simulation results

% Tha path to save the simulation results to
results_path = simulationDir(directory_label);

if ~exist(results_path, 'dir')
    mkdir(results_path)
end

% Are we running in GNU Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% Add the required libraries if we are running in octave
if isOctave
    pkg load statistics;
    pkg load image;
end

addpath(results_path);

%% Sample import and plotting
% Importing the sample as a TriagSurface object.
% Third argument asks if we wish to plot the surface in 3D, will save it to the
% simultion path.
switch sample_type
    case 'flat'
        sample_surface = flatSample(square_size, dist_to_sample, defMaterial);
        make_sphere = 0;
        sample_description = ['A flat square sample size ' ...
            num2str(square_size) 'mm.'];
    case 'strips'
        sample_surface = strip_series(dist_to_sample, working_dist);
        make_sphere = 0;
        sample_description = 'A sample made up of two series of strips with properties varying across';
    case 'sphere'
        sample_surface = flatSample(square_size, dist_to_sample, defMaterial);
        make_sphere = 1;
        sample_description = ['A single analytic sphere, radius ' ...
            num2str(sphere_r) 'mm on a flat square of ' num2str(square_size) 'mm.'];
    case 'custom'
        sample_surface = inputSample('fname', sample_fname, 'sampleDist', dist_to_sample, ...
                                     'workingDist', working_dist, 'scale', scale, ...
                                     'defMaterial', defMaterial);
        make_sphere = 0;
    case 'photoStereo'
        sample_surface = photo_stereo_test(working_dist);
        make_sphere = 1;
        sphere_r = 0.05;
        sphere_c = [-0.1, -dist_to_sample - sphere_r*2/3, -0.1];
        diffuse = [1, 90*pi/180];
    case 'special'
        sample_surface = inputSample('fname', sample_fname, 'dontMeddle', true, 'scale', 10e-4);
        sample_surface.rotateZ;
        sample_surface.moveBy([0, -2.121, 0]);
        make_sphere = 0;
end

sample_surface.reflect_axis('x');

if strcmp(typeScan, 'line')
    sample_surface.moveBy(init_displacement);
end

% A struct to represent the sphere
sphere_c = [dist_to_sample*tand(init_angle), dist_to_sample + sphere_r, 0];
sphere.c = sphere_c;
sphere.make = make_sphere;
sphere.r = sphere_r;
sphere.material = defMaterial;

% Do any extra manipulation of the sample here
% if false
%     %sample_surface.rotateX;
%     sample_surface.rotateY;
%     sample_surface.rotateY;
%     %sample_surface.rotateY;
%     %sample_surface.moveBy([0, -2.2, 2.3]);
% end

% Plot the sample surface in 3D space, if we are using a graphical window
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
        print(fullfile(results_path, 'sample_closeUp.eps'), '-depsc');
    end
end

%% Pinhole plate import and plotting
switch pinhole_model
    case 'stl'
        pinhole_surface = import_plate(plate_accuracy);

        % Plot if using a graphical window
        if ~isOctave
            if feature('ShowFigureWindows')
                sample_surface.patchPlot(true);
                pinhole_surface.patchPlot(false);
                view([-5 -5 5]);
            end
        else
            sample_surface.patchPlot(true);
            pinhole_surface.patchPlot(false);
            view([-5 -5 5]);
        end

        % To pass to the functions
        thePlate = 0;
        apertureAbstract = 0;
    case 'new'
        pinhole_surface = import_newPlate(plate_accuracy);

        % Plot if using a graphical window
        if ~isOctave
            if feature('ShowFigureWindows')
                sample_surface.patchPlot(true);
                pinhole_surface.patchPlot(false);
                view([-5 -5 5])
            end
        else
            sample_surface.patchPlot(true);
            pinhole_surface.patchPlot(false);
            view([-5 -5 5])
        end

        % To pass to the functions
        thePlate = 0;
        apertureAbstract = 0;

        pinhole_model = 'stl';
    case 'new_micro'
        % TODO
        error('Not written this bit of code yet...');
    otherwise
        % Do not have a CAD repreentation of the pinhole plate

        % Create an empty TriagSurface as the pinhole plate
        pinhole_surface = TriagSurface();

        % List with the information about the plate in
        % TODO: use a struct rather than a cell array.
        thePlate = {plate_represent, n_detectors, circle_plate_r, aperture_axes, aperture_c};
        apertureAbstract = {aperture_theta, aperture_phi, aperture_half_cone};
end

%% Compile the mex files

if recompile
    mexCompile();
end

% traceSimpleMultiGen('sample', sample_surface, 'maxScatter', maxScatter, 'plate', thePlate, 'dist', dist_to_sample, 'sphere', sphere, 'source', direct_beam.source_model, 'beam', direct_beam);
% return;

%% Performing the simulation
switch typeScan
    case 'rectangular'
        % For a rectangular scan
        simulationData = rectangularScan(sample_surface, xrange, zrange, ...
            direct_beam, raster_movment2D_x, raster_movment2D_z, ...
            maxScatter, pinhole_surface, effuse_beam, ...
            dist_to_sample, sphere, results_path, pinhole_model, ...
            thePlate, apertureAbstract, ray_model, n_detectors);
    case 'multiple_rectangular'
        simulationData = multipleRectangularScan(sample_surface, range_y, raster_movement_y,...
            xrange, zrange, direct_beam, raster_movment2D_x, raster_movment2D_z, ...
            maxScatter, pinhole_surface, effuse_beam, ...
            dist_to_sample, sphere, results_path, pinhole_model, ...
            thePlate, apertureAbstract, ray_model, n_detectors);
    case 'line'
        % For a line scan
        % TODO: update with the new lower level functions
        simulationData = lineScan(sample_surface, range1D, direct_beam, ...
            raster_movment1D, maxScatter, Direction, pinhole_surface, effuse_beam, ...
            dist_to_sample, sphere, results_path, pinhole_model, ...
            thePlate, apertureAbstract, ray_model);
    case 'single pixel'
        % For a single pixel
        % TODO: update with the new lower level functions
        simulationData = singlePixel(sample_surface, direct_beam, ...
            maxScatter, pinhole_surface, effuse_beam, dist_to_sample, sphere, ...
            results_path, save_to_text, pinhole_model, thePlate, apertureAbstract);
    case 'rotations'
        % Perform multiple scans while rotating the sample in between.
        simulationData = {};
        h = waitbar(0, 'Proportion of simulations performed');
        N = length(rot_angles);
        sphere_centre = sphere.c;

        % Loop through the rotations
        for i_=1:N
            s_surface = copy(sample_surface);
            s_surface.rotateGeneral('y', rot_angles(i_));

            % Rotate the centre of the sphere
            theta = rot_angles(i_)*pi/180;
            s = sin(theta);
            c = cos(theta);
            R = [c, 0, s; 0, 1, 0; -s, 0, c];
            sphere.c = (R*sphere_centre')';

            subPath = [results_path '/rotation' num2str(rot_angles(i_))];
            if ~exist(subPath, 'dir')
                mkdir(subPath)
            end

            simulationData{i_} = rectangularScan(s_surface, xrange, zrange, ...
                direct_beam, raster_movment2D_x, raster_movment2D_z, ...
                maxScatter, pinhole_surface, effuse_beam, ...
                dist_to_sample, sphere, subPath, pinhole_model, ...
                thePlate, apertureAbstract, ray_model, n_detectors); %#ok<SAGROW>

            waitbar(i_/N, h);
        end
        close(h);
        delete(h);
    otherwise
        error(['Need to specify a valid type of scan: "line", ', ...
               '"rectangular", "single pixel"']);
end

% kill parallel pool
% delete(gcp('nocreate'));

%% Output data about simulation to files

% Create input structs to hole input data

scan_inputs.type_scan = typeScan;
scan_inputs.maxScatter = maxScatter;
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
        scan_inputs.raster_movment1D = raster_movment1D;
        scan_inputs.range1D = range1D;
        scan_inputs.direction_1D = Direction;
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
sample_inputs.material = defMaterial;
sample_inputs.dist_to_sample = dist_to_sample;
sample_inputs.sphere = sphere;
sample_inputs.sample_description = sample_description;
if strcmp(sample_type, 'rectangular')
    sample_inputs.square_size = square_size;
else
    sample_inputs.square_size = NaN;
end

pinhole_plate_inputs.pinhole_model = pinhole_model;
pinhole_plate_inputs.working_dist = working_dist;
switch pinhole_model
    case {'stl', 'new'}
        pinhole_plate_inputs.plate_accuracy = plate_accuracy;
        pinhole_plate_inputs.n_detectors = 1;
        pinhole_plate_inputs.plate_represent = 1;
        pinhole_plate_inputs.aperture_axes = NaN;
        pinhole_plate_inputs.aperture_c = NaN;
    case 'circle'
        pinhole_plate_inputs.plate_accuracy = NaN;
        pinhole_plate_inputs.n_detectors = 1;
        pinhole_plate_inputs.plate_represent = plate_represent;
        pinhole_plate_inputs.aperture_axes = aperture_axes;
        pinhole_plate_inputs.aperture_c = aperture_c;
    case 'N circle'
        pinhole_plate_inputs.plate_accuracy = NaN;
        pinhole_plate_inputs.n_detectors = n_detectors;
        pinhole_plate_inputs.plate_represent = plate_represent;
        pinhole_plate_inputs.aperture_axes = aperture_axes;
        pinhole_plate_inputs.aperture_c = aperture_c;
    case 'abstract'
        % TODO: make the 'abstract' simulations work
        pinhole_plate_inputs.n_detectors = n_detectors;
        pinhole_plate_inputs.plate_represent = 0;
end

% Save all data to a .mat file
if output_data
    save(fullfile(results_path, data_fname), 'simulationData', 'sample_inputs', ...
        'direct_beam', 'effuse_beam', 'pinhole_plate_inputs', 'scan_inputs');
end

% Save formatted data to a .mat file, does not include all the parameters
% but does include the core outputs.
if output_data && strcmp(typeScan, 'rotations')
    simulationData.formatOutput(rot_angles, results_path);
elseif output_data
    simulationData.formatOutput(results_path);
end

% Save main data to a text file
textFname = 'data_for_plotting.csv';
if save_to_text && strcmp(typeScan, 'rotations')
    for i_=1:length(rot_angles)
        currentFname = [textFnam(1:end-4) num2str(rot_angles(i_)) '.csv'];
        simulationData.saveText([results_path '/' currentFname]);
    end
elseif save_to_text
    simulationData.saveText([results_path '/' textFname]);
end

% Save the parameters to a text file
if saveParams
    saveParam(sample_inputs, direct_beam, effuse_beam, ...
        pinhole_plate_inputs, scan_inputs);
end

