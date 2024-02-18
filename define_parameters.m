% Reads parameters from the parameter text file and sorts them out into
% useful data structures to pass to the main simulation. There are 2 basic
% sections, the first reads parameters from the text file, the second turns
% defines other parameters not in the text file, and the third puts the
% parameters into input structs.

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
aperture_rotate = parse_list_input(param_list{7});
rot_angles = parse_rotations(param_list{8});
det_pos = parse_rotations(param_list{9});
pinhole_model = parse_pinhole(param_list{10});
aperture_size = parse_list_input(param_list{11});
aperture_theta = aperture_size(1);
aperture_phi = aperture_size(2);
aperture_half_cone = str2double(param_list{12});

% Set up source
n_rays = str2double(param_list{13});
pinhole_r = str2double(param_list{14});
source_model = strtrim(param_list{15});
theta_max = str2double(param_list{16});
sigma_source = str2double(param_list{17});
if ~parse_yes_no(param_list{18})
    effuse_size = 0;
else
    effuse_size = str2double(param_list{19});
end

% Set up sample
sample_type = strtrim(param_list{20});
material = parse_scattering(strtrim(param_list{21}), str2double(param_list{22}), ...
    str2double(param_list{23}));
sample_description = param_list{24};
dist_to_sample = str2double(param_list{25});
sphere_rs = parse_list_input(param_list{26});
n_sphere = length(sphere_rs);
sphere_cs = reshape(parse_list_input(param_list{27}), [2, n_sphere]);
circle_r = str2double(param_list{28});
square_size = str2double(param_list{29});
sample_fname = strtrim(param_list{30});
dontMeddle = parse_yes_no(param_list{31});

% Set up scan
pixel_seperation = str2double(param_list{32});
range_x = str2double(param_list{33});
range_z = str2double(param_list{34});
init_angle_pattern = ~parse_yes_no(param_list{35});

% 1D scan parameters
scan_direction = strtrim(param_list{36});
range_1D = [str2double(param_list{37}), str2double(param_list{38})];
res_1D = str2double(param_list{39});

% Other parameters
directory_label = strtrim(param_list{40});
recompile = parse_yes_no(param_list{41});

%% Generate parameters from the inputs 

pinhole_c = [-working_dist*tand(init_angle), 0, 0];
n_effuse = n_rays*effuse_size;
raster_movment2D_x = pixel_seperation;
raster_movment2D_z = pixel_seperation;
xrange = [-range_x/2, range_x/2];% + tand(init_angle)*(dist_to_sample - working_dist);
zrange = [-range_z/2, range_z/2];

% TODO: sphere locations and 
sphere_y = -dist_to_sample + sphere_rs;
sphere_cs = [sphere_cs(1,:); sphere_y; sphere_cs(2,:)];

%% Define remaining parameters

circle_c = [0, -working_dist, 0];
circle_n = [0, 1, 0];

% The maximum number of sample scatters per ray. There is a hard-coded total
% maximum number of scattering events of 1000 (sample and pinhole plate). Making
% this uneccaserily large will increase the memory requirments of the
% simulation.
max_scatter = 100;

% If rotations are present the scan pattern can be regular or be adjusted to
% match the rotation of the sample
scan_pattern = 'rotations';

% Do we want to generate rays in Matlab (more flexibility, more output options)
% or in C (much lower memory requirments and slightly faster), 'C' or 'MATLAB'
% In general stick to 'C' unless your own source model is being used
ray_model = 'C';
if strcmp(source_model, 'Infinite')
    ray_model = 'C';
    pinhole_r = 0;
end

% Exponant of the cosine in the effuse beam model
cosine_n = 1;

% In the case of the predefined CAD model, specify the accuraccy of the
% triangulation, 'low', 'medium', or 'high' (always use 'low').
plate_accuracy = 'low';

% In the case of 'circle', specify the radius of the circle (mm).
circle_plate_r = 4;

% Should a flat pinhole plate be modelled (with 'N circle'). not including may
% speed up the simulation but won't model the effuse and multiple scattering
% backgrounds properly.
plate_represent = 0;


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
scale = 0.5;
% By default have the scale set to 1
scale = 1;
%scale = 2;

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

%% Create input structs to hold input data

% TODO: convert from structs into classes?? <- generally everything would
% be better than way, but ideally most of the code would be migrated to
% C...
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
    case 'line_detector_slide'
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
        pinhole_plate_inputs.aperture_rotate = NaN;
        pinhole_plate_inputs.circle_plate_r = NaN;
        pinhole_plate_inputs.aperture_theta = NaN;
        pinhole_plate_inputs.aperture_phi = NaN;
        pinhole_plate_inputs.aperture_half_cone = NaN;
    case 'N circle'
        pinhole_plate_inputs.plate_accuracy = NaN;
        pinhole_plate_inputs.n_detectors = n_detectors;
        pinhole_plate_inputs.plate_represent = plate_represent;
        pinhole_plate_inputs.aperture_axes = aperture_axes;
        pinhole_plate_inputs.aperture_c = aperture_c;
        pinhole_plate_inputs.aperture_rotate = aperture_rotate;
        pinhole_plate_inputs.circle_plate_r = circle_plate_r;
        pinhole_plate_inputs.aperture_theta = NaN;
        pinhole_plate_inputs.aperture_phi = NaN;
        pinhole_plate_inputs.aperture_half_cone = NaN;
    case 'abstract'
        pinhole_plate_inputs.plate_accuracy = NaN;
        pinhole_plate_inputs.n_detectors = 1;
        pinhole_plate_inputs.plate_represent = 0;
        pinhole_plate_inputs.aperture_axes = NaN;
        pinhole_plate_inputs.aperture_c = NaN;
        pinhole_plate_inputs.aperture_rotate = NaN;
        pinhole_plate_inputs.circle_plate_r = NaN;
        pinhole_plate_inputs.aperture_theta = aperture_theta;
        pinhole_plate_inputs.aperture_phi = aperture_phi;
        pinhole_plate_inputs.aperture_half_cone = aperture_half_cone;
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