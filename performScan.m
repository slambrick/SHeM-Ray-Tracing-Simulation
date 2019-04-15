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
maxScatter = 10;

% Type of scan 'line', 'rectangular', or 'single pixel'
typeScan = 'rectangular';

% Recompile mex files?
% Required if using on a new computer or if changes to .c files have been made.
recompile = true;

%% Beam/source parameters %%

% The inicidence angle in degrees
init_angle = 45;

% Geometry of pinhole
pinhole_c = [-2.121, 0, 0];
pinhole_r = 0.0025;

% Number of rays to use and the width of the source
n_rays = 10000;

% skimmer radius over source - pinhole distance
theta_max = atan(0.05/100); 

% Standard deviationo of the Gaussian model of the source
sigma_source = 0.1;

% Model for the virtual source to use, 'Gaussian' or 'Effuse'
source_model = 'Uniform';

% Put the information in a cell array to pass through functions
% TODO: use a struct rather than a cell array.
direct_beam = {n_rays, pinhole_c, pinhole_r, theta_max, source_model, ...
    init_angle, sigma_source};

% Do we want to generate rays in Matlab (more flexibility, more output options)
% or in C (much lower memory requirments and slightly faster), 'C' or 'MATLAB'
ray_model = 'C';

%% Effusive beam parameters

% If a cosine is being used then specify the exponent of the cosine
cosine_n = 1;

% How large is the effusive beam (proportion of the size of the main beam). Set
% to zero if the effusive beam is not to be moddeled
effuse_size = 10;

% Information on the effuse beam
n_effuse = n_rays*effuse_size;
% TODO: use a struct rather than a cell array.
effuse_beam = {n_effuse, pinhole_c, pinhole_r, cosine_n};

%% Pinhole plate parameters
% Specify how to model the pinhole plate:
%  'stl'       - Use the predefined CAD model of the pinhole plate (plate as it
%                is Feb 2018)
%  'aperture'  - Use only the detector aperture, fastest but ignores the multiple
%                scattering and will under-represent the effuse contribution.
%                Both this option and 'circle' will not included transmission
%                proabability through the detector cone
%  'circle'    - Use the detector aperture and model the pinhole plate as a
%                circle
%  'N circle'  -
%  'new'       - TODO
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
% direction ('z'). The apert0ure is always centred on the x-axis and is displaced
% by the specified amount.
n_detectors = 3;
aperture_axes = [1, 1, 1, 1, 1, 1];
aperture_c = [2, 0, 1, 2, 1, -2];
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
raster_movment2D_x = 0.02*sqrt(2);
raster_movment2D_z = 0.02;
xrange = [-0.3, 0.3];
zrange = [-0.2, 0.2];

%% Parameters for a 1d scan
% For line scans in the y-direction be careful that the sample doesn't go
% behind the pinhole plate.
raster_movment1D = 1e-2;
range1D = [-1 3];
Direction = 'y';

%% Sample parameters
% The sample file, include the full path
sample_fname = 'simulations/001_example/example_sample.stl';

% Sample scaling, for if the CAD model had to be made at a larger scale. 10 will
% make the model 10 times larger (Inventor exports in cm by default...).
scale = 0.2;

% A string giving a brief description of the sample, for use with
% sample_type = 'custom'
sample_description = 'A sample with some example topography.';

% What type of sample to use :
%  'flat'   - A flat square (need to specify square_size)
%  'sphere' - An anlaytic sphere on a flat square surface (need to specify 
%             square_size and sphere_r)
%  'custom' - Uses the CAD model provided in the .stl file
sample_type = 'sphere';

% The level of diffuse scattering for the sample, between 0 and 1, 2 gives a 
% uniform distribution. This is used for both triangulated surface and the
% analytic sphere if it is being used.
diffuse = [1, 90*pi/180];

% How close should the nearest point of the sample be to the pinhole plate, the
% defualt is 2.121 to maintain the 45o geometry. If an analytic sphere is being
% used then this is the distance between the flat surface the sphere sits on and
% the pinhole plate.
dist_to_sample = 2;

% The radius of the anayltic sphere (mm) (if it being included)
sphere_r = 0.05;

% If a flat sample is being used or if a sphere on a flat surface what is the 
% length of the sides of the square.
square_size = 4;

%% Output and plotting parameters

% Where to save figures/data files
% All figures and output data will be saved to this directory.
directory_label = 'doesitwork';

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
        error('Specify a correct model of pinhole plate.');
    end
end

%% Paths to functions
addpath('stlread', 'functions', 'functions/interface_functions', 'classes', ...
        'mexFiles', 'DylanMuir-ParforProgMon-a80c9e9', 'functions/standard_samples');

%% Path for simulation results

% Tha path to save the simulation results to
thePath = simulationDir(directory_label);

if ~exist(thePath, 'dir')
    mkdir(thePath)
end

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
            'scattering_parameters', diffuse(2));
        make_sphere = 0;
end

% TODO: use a struct rather than a cell array.
sphere = {make_sphere, sphere_r, diffuse(1), diffuse(2)};

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
    
    if strcmp(typeScan, 'rectangular') 
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
        pinhole_surface = import_newPlate();

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
% must include the GSL libraries
% The -ffast-math improves speed
if recompile
    % For stl pinhole plate and sample
    mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
        mexFiles/tracingMex.c            mexFiles/trace_ray.c ...
        mexFiles/tracing_functions.c ...
        mexFiles/scattering3D.c          mexFiles/scattering_processes3D.c ...
        mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
        mexFiles/common_helpers.c
    mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
        mexFiles/tracingGenMex.c         mexFiles/trace_ray.c ...
        mexFiles/tracing_functions.c ...
        mexFiles/scattering3D.c          mexFiles/scattering_processes3D.c ...
        mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
        mexFiles/common_helpers.c
    
    % For a simple model of the pinhole plate
    mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
        mexFiles/tracingSimpleMex.c      mexFiles/trace_ray.c ...
        mexFiles/tracing_functions.c ...
        mexFiles/scattering3D.c          mexFiles/scattering_processes3D.c ...
        mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
        mexFiles/common_helpers.c
    mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
        mexFiles/tracingSimpleGenMex.c   mexFiles/trace_ray.c ...
        mexFiles/tracing_functions.c ...
        mexFiles/scattering3D.c          mexFiles/scattering_processes3D.c ...
        mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
        mexFiles/common_helpers.c
    mex -lgsl -lm -lgslcblas CFLAGS="\$CFLAGS -Wall" ...
        mexFiles/tracingMultiGenMex.c   mexFiles/trace_ray.c ...
        mexFiles/tracing_functions.c ...
        mexFiles/scattering3D.c          mexFiles/scattering_processes3D.c ...
        mexFiles/ray_tracing_structs3D.c mexFiles/small_functions3D.c ...
        mexFiles/common_helpers.c

    %mex -lgsl -lgslcblas x mexFiles/noPlateMex.c ...
    %    mexFiles/tracing_functions.c mexFiles/small_functions.c
    %mex -lgsl -lgslcblas CFLAGS="\$CFLAGS -ffast-math" mexFiles/cosineMex.c ...
    %    mexFiles/small_functions.c
    %mex -lgsl -lgslcblas CFLAGS="\$CFLAGS -ffast-math" mexFiles/abstractPlateMex.c ...
    %    mexFiles/tracing_functions.c mexFiles/small_functions.c
    
    % A binning function I wrote.
    mex mexFiles/binMyWayMex.c
end

%% Performing the simulation
switch typeScan
    case 'rectangular'
        % For a rectangular scan
        simulationData = rectangularScan(sample_surface, xrange, zrange, ...
            direct_beam, raster_movment2D_x, raster_movment2D_z, ...
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
    otherwise
        error(['Need to specify a valid type of scan: "line", ', ...
               '"rectangular", "single pixel"']);
end

%% Output data about simulation to files
if output_data
    save([thePath '/' data_fname]);
end

% Save the parameters to a text file
% TODO: do this better
if false %saveParams 
    fid = fopen([thePath '/' paramsFile], 'w');
    
    FORMAT1 = '%s = ';
    FORMAT2 = '%2.8f\n';
    FORMAT3 = '%i\n';
    
    fprintf(fid, '%s\n\n', 'Parameters and results from simulation');
    
    fprintf(fid, '%s %s %s\n\n', 'This is a', typeScan, 'scan.');
    
    fprintf(fid, '%s\n\n', sample_description);
    
    if strcmp(sample_type, 'custom')
        fprintf(fid, FORMAT1, 'The scaling used for the sample model');
        fprintf(fid, FORMAT2, scale);
    end
    
    fprintf(fid, FORMAT1, 'Date');
    fprintf(fid, '%s\n', datestr(now, 'dd-mm-yyyy'));
    
    fprintf(fid, FORMAT1, 'Minimum distance from the pinholePlate to the sample (mm)');
    fprintf(fid, FORMAT2, dist_to_sample);
    
    fprintf(fid, FORMAT1, 'Scattering label for the sample');
    fprintf(fid, FORMAT2, diffuse);
    
    fprintf(fid, FORMAT1, 'Ray seperation (mm)');
    fprintf(fid, FORMAT2, ray_sep);
    
    fprintf(fid, FORMAT1, 'Maximum number of allowed Scatters');
    fprintf(fid, FORMAT3, maxScatter);
    
    fprintf(fid, FORMAT1, 'Pinhole radius (mm)');
    fprintf(fid, FORMAT2, pinhole_r);
    
    fprintf(fid, FORMAT1, 'Pinhole center (x, z) (mm)');
    fprintf(fid, '%f, %f\n', pinhole_c(1), pinhole_c(3));
    
    fprintf(fid, FORMAT1, 'Total number of rays per pixel');
    fprintf(fid, FORMAT3, n_rays);
    
    fprintf(fid, FORMAT1, 'Number of rays per point in the pinhole');
    fprintf(fid, FORMAT3, multipl);
    
    fprintf(fid, FORMAT1, 'Relative size of the effuse beam');
    fprintf(fid, FORMAT3, effuse_size);
    
    fprintf(fid, FORMAT1, 'Model used for the pinhole plate');
    fprintf(fid, '%s\n', pinhole_model);
    
    switch pinhole_model
        case 'stl'
            fprintf(fid, FORMAT1, 'Level of accuracy of pinhole plate');
            fprintf(fid, '%s\n', plate_accuracy);
        case 'aperture'
            fprintf(fid, FORMAT1, 'Centre of the detector aperture in x (mm)');
            fprintf(fid, FORMAT2, aperture_c);
            fprintf(fid, FORMAT1, 'Axes of the aperture (x,z/mm)');
            fprintf(fid, '%f, %f\n', aperture_axis_1, aperture_axis_2);
        case 'circle'
            fprintf(fid, FORMAT1, 'Centre of the detector aperture in x (mm)');
            fprintf(fid, FORMAT2, aperture_c);
            fprintf(fid, FORMAT1, 'Axes of the aperture (x,z/mm)');
            fprintf(fid, '%f, %f\n', aperture_axis_1, aperture_axis_2);
            fprintf(fid, FORMAT1, 'Diameter of the circular pinhole plate (mm)');
            fprintf(fid, FORMAT2, circle_plate_r);
    end
    
    if strcmp(typeScan, 'rectangular')
        fprintf(fid, FORMAT1, 'Number of pixels in x');
        fprintf(fid, FORMAT3, simulationData.nx_pixels);
        
        fprintf(fid, FORMAT1, 'Number of pixels in z');
        fprintf(fid, FORMAT3, simulationData.nz_pixels);
        
        fprintf(fid, FORMAT1, 'Seperation of pixels (mm)');
        fprintf(fid, FORMAT2, raster_movment2D_z);
        
        fprintf(fid, FORMAT1, 'Range of x values (min, max) (mm)');
        fprintf(fid, '%f, %f\n', xrange(1), xrange(2));
        
        fprintf(fid, FORMAT1, 'Range of z values (min, max) (mm)');
        fprintf(fid, '%f, %f\n', zrange(1), zrange(2));
    end
    
    if strcmp(typeScan, 'line')
        fprintf(fid, FORMAT1, 'Seperation of pixels (mm)');
        fprintf(fid, FORMAT2, raster_movment1D);
        
        fprintf(fid, '%s = %s\n', 'Direction of line scan', Direction);
        
        fprintf(fid, FORMAT1, 'Range of positions (min, max) (mm)');
        fprintf(fid, '%f, %f\n', range1D(1), range1D(2));
    end
        
    if ~strcmp(typeScan, 'single pixel')
        fprintf(fid, FORMAT1, 'Total number of pixels');
        fprintf(fid, FORMAT3, simulationData.N_pixels);
        
        fprintf(fid, FORMAT1, 'Estimate for the time taken (s)');
        fprintf(fid, FORMAT2, simulationData.time_estimate);
    end
    
    fprintf(fid, FORMAT1, 'Time the core simulation took (s)');
    fprintf(fid, FORMAT2, simulationData.time);
    
    fclose(fid);
end

