% Copyright (c) 2018-20, Sam Lambrick.
% All rights reserved.
% This file is part of the SHeM Ray Tracing Simulation, subject to the 
% GNU/GPL-3.0-or-later.

clear
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Paths to functions
addpath('stlread', genpath('functions'), 'classes', ...
        'mexFiles', 'DylanMuir-ParforProgMon-a80c9e9', ...
        'surf2stl');

%% Read parameters from text file
param_fname = 'ray_tracing_parameters.txt';
param_list = read_parameters(param_fname);

% Set up virtual microscope
working_dist = str2double(param_list{1});
init_angle = str2double(param_list{2});
typeScan = strtrim(param_list{3});
n_detectors = str2double(param_list{4});
aperture_axes = parse_list_input(param_list{5});
aperture_c = parse_list_input(param_list{6});
rot_angles = parse_list_input(param_list{7});

% Set up source
n_rays = str2double(param_list{8});
pinhole_r = str2double(param_list{9});
source_model = strtrim(param_list{10});
theta_max = str2double(param_list{11});
sigma_source = str2double(param_list{12});
if ~parse_yes_no(param_list{13})
    effuse_size = 0;
else
    effuse_size = str2double(param_list{14});
end

% Set up sample
sample_type = strtrim(param_list{15});
diffuse = parse_scattering(strtrim(param_list{16}), str2double(param_list{17}), ...
    str2double(param_list{18}));
sample_description = param_list{19};
dist_to_sample = str2double(param_list{20});
sphere_r = str2double(param_list{21});
square_size = str2double(param_list{22});
sample_fname = strtrim(param_list{23});

% Set up scan
pixel_seperation = str2double(param_list{24});
range_x = str2double(param_list{25});
range_z = str2double(param_list{26});
init_angle_pattern = ~parse_yes_no(param_list{27});

% Other parameters
directory_label = strtrim(param_list{28});
recompile = parse_yes_no(param_list{29});

%% Generate parameters from the inputs 

pinhole_c = [-working_dist*tand(init_angle), 0, 0];
n_effuse = n_rays*effuse_size;
raster_movment2D_x = pixel_seperation;
raster_movment2D_z = pixel_seperation;
xrange = [-range_x/2, range_x/2];
zrange = [-range_z/2, range_z/2];
sphere_c = [0, -dist_to_sample, 0];

%% Define remaining parameters

% The maximum number of sample scatters per ray. There is a hard-coded total
% maximum number of scattering events of 1000 (sample and pinhole plate). Making
% this uneccaserily large will increase the memory requirments of the
% simulation.
maxScatter = 20;

% If rotations are present the scan pattern can be regular or be adjusted to
% match the rotation of the sample
scan_pattern = 'regular';

% Do we want to generate rays in Matlab (more flexibility, more output options)
% or in C (much lower memory requirments and slightly faster), 'C' or 'MATLAB'
ray_model = 'C';

% Exponant of the cosine in the effuse beam model
cosine_n = 1;

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

% Should a flat pinhole plate be modelled (with 'N circle'). not including may
% speed up the simulation but won't model the effuse and multiple scattering
% backgrounds properly.
plate_represent = 1;

% In the case of 'abstract', specify the two angles of the location of the
% detector aperture and the half cone angle of its extent. Note that the
% aperture can only be placed in the hemisphere facing the sample. All
% angles in degrees.
% TODO: this
aperture_theta = 0;
aperture_phi = 0;
aperture_half_cone = 15;

% For line scans in the y-direction be careful that the sample doesn't go
% behind the pinhole plate.
raster_movment1D = 200e-3;
range1D = [-1 4];
Direction = 'y';

% Sample scaling, for if the CAD model had to be made at a larger scale. 10 will
% make the model 10 times smaller (Inventor exports in cm by default...).
scale = 2;

%% Create parameter structs

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

%% Output and plotting parameters

% Which figures to plot
% The starting positions of the rays and the number of rays at each point
plot_density = false;

% Save simulation data to file
output_data = true;
data_fname = 'scatteringData.mat';
saveParams = true;
paramsFile = 'scatteringParameters.txt';
% Applies only to line scans, saves in thePath with name 'data_for_plotting.csv'
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
if isOctave
    pkg load statistics;
    pkg load image;
end

addpath(thePath);

%% Sample import and plotting
% Importing the sample as a TriagSurface object.
% Third argument asks if we wish to plot the surface in 3D, will save it to the
% simultion path.
switch sample_type
    case 'flat'
        sample_surface = flatSample(square_size, dist_to_sample, diffuse(1));
        make_sphere = 0;
        sample_description = ['A flat square sample size ' ...
            num2str(square_size) 'mm.'];
    case 'sphere'
        sample_surface = flatSample(square_size, dist_to_sample, diffuse(1), diffuse(2));
        make_sphere = 1;
        sample_description = ['A single analytic sphere, radius ' ...
            num2str(sphere_r) 'mm on a flat square of ' num2str(square_size) 'mm.'];
    case 'custom'
        sample_surface = inputSample('fname', sample_fname, 'scattering', ...
            diffuse(1), 'plate_dist', dist_to_sample, 'scale', scale, ...
            'parameters', diffuse(2), 'working_dist', working_dist);
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

% A struct to represent the sphere
shere_c = [dist_to_sample*tand(init_angle), dist_to_sample + sphere_r, 0];
sphere.make = make_sphere;
sphere.r = sphere_r;
sphere.scattering = diffuse(1);
sphere.scattering_parameter = diffuse(2);
sphere.c = sphere_c;

% Do any extra manipulation of the sample here
if false
    %sample_surface.rotateX;
    sample_surface.rotateY;
    sample_surface.rotateY;
    %sample_surface.rotateY;
    %sample_surface.moveBy([0, -2.2, 2.3]);
end

% Plot the sample surface in 3D space, if we are using a graphical window
if true %feature('ShowFigureWindows')
    if ~strcmp(typeScan, 'single_pixel')
        sample_surface.patchPlot(true);
        ylim([-dist_to_sample - 0.2, -dist_to_sample + 0.2]);
    end
    
    if strcmp(typeScan, 'rectangular') || strcmp(typeScan, 'rotations')
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
        print([thePath '/sample_closeUp.eps'], '-depsc');
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
        thePlate.plate_represent = plate_represent;
        thePlate.n_detectors = n_detectors;
        thePlate.circle_plate_r = circle_plate_r;
        thePlate.aperture_axes = aperture_axes;
        thePlate.aperture_c = aperture_c;
        %thePlate = {plate_represent, n_detectors, circle_plate_r, aperture_axes, aperture_c};
        apertureAbstract = {aperture_theta, aperture_phi, aperture_half_cone};
end
    
%% Compile the mex files

files_exist = exist('tracingMultiGenMex.mexa64', 'file');

if recompile || ~files_exist
    mexCompile();
end

%% Create input structs to hole input data

scan_inputs.type_scan = typeScan;
scan_inputs.maxScatter = maxScatter;
switch typeScan
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
        scan_inputs.direction_1D = direction;
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
sample_inputs.scattering = diffuse;
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


%% Performing the simulation
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

        simulationData = rectangularScan(sample_surface, raster_pattern, ...
            direct_beam, ...
            maxScatter, pinhole_surface, effuse_beam, ...
            dist_to_sample, sphere, thePath, pinhole_model, ...
            thePlate, apertureAbstract, ray_model, n_detectors);
    case 'line'
        % For a line scan
        % TODO: update with the new lower level functions
        simulationData = lineScan(sample_surface, range1D, direct_beam, ...
            raster_movment1D, maxScatter, Direction, pinhole_surface, effuse_beam, ...
            dist_to_sample, sphere, thePath, save_to_text, pinhole_model, ...
            thePlate, apertureAbstract, ray_model);
    case 'single pixel'
        % For a single pixel
        % TODO: update with the new lower level functions
        simulationData = singlePixel(sample_surface, direct_beam, ...
            maxScatter, pinhole_surface, effuse_beam, dist_to_sample, sphere, ...
            thePath, save_to_text, pinhole_model, thePlate, apertureAbstract);
    case 'rotations'
        % Perform multiple scans while rotating the sample in between.
        simulationData = {};
        h = waitbar(0, 'Proportion of simulations performed', 'Name', 'Ray tracing progress');
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
            
            simulationData{i_} = rectangularScan(s_surface, raster_pattern, ...
                direct_beam, ...
                maxScatter, pinhole_surface, effuse_beam, ...
                dist_to_sample, sphere, subPath, pinhole_model, ...
                thePlate, apertureAbstract, ray_model, n_detectors); %#ok<SAGROW>
            
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
    otherwise
        error(['Need to specify a valid type of scan: "line", ', ...
               '"rectangular", "single pixel"']);
end

%% Output data about simulation to files

% Save all data to a .mat file
if output_data
    save([thePath '/' data_fname], 'simulationData', 'sample_inputs', ...
        'direct_beam', 'effuse_beam', 'pinhole_plate_inputs', 'scan_inputs');
end

% Save formatted data to a .mat file, does not include all the parameters
% but does include the core outputs.
if output_data && strcmp(typeScan, 'rotations')
    formatOutputRotation(simulationData, thePath, rot_angles);
elseif output_data
    simulationData.formatOutput(thePath);
end

% Save main data to a text file
textFname = 'data_for_plotting.csv';    
if save_to_text && strcmp(typeScan, 'rotations')
    for i_=1:length(rot_angles)
        % TODO: bug here
        currentFname = [textFname(1:end-4) num2str(rot_angles(i_)) '.csv'];
        simulationData.saveText([thePath '/' currentFname]);
    end
elseif save_to_text
    simulationData.saveText([thePath '/' textFname]);
end

% Save the parameters to a text file
%if saveParams
%    saveParam(sample_inputs, direct_beam, effuse_beam, ...
%        pinhole_plate_inputs, scan_inputs);
%end

