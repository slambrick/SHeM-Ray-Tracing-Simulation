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
maxScatter = 100;

% Type of scan 'line', 'rectangular', or 'single pixel'
typeScan = 'rectangular';

% Recompile mex files?
% Required if using on a new computer or if changes to .c files have been made.
recompile = true;

%% Beam/source parameters

% Geometry of pinhole. Specify the centre of the pinhole and the maximum radius
% of the beam that is to be created at the pinhole. When using a directional
% distribution to model the spread of the beam the radius here is the true
% radius of the pinhole, when using a parallel beam with a density profile it is
% the maximum extent of the  
pinhole_c = [-2.121, 0, 0];
pinhole_r = 0.005;

% Max number of rays at each location.
% 30 produces decent 2D images.
% Making this number too high will make the simulation slow and increase the
% memory requirements (current configuration 40-80 should work).
multipl = 40;

% Seperation of the rays (mm), they are arranged in a square grid across the
% whole of the pinhole beam.
% 0.0001mm is reasonable for an image.
ray_sep = 0.0001;

% Should the rays be generated all parallel (true) with a density profile that
% is representative of the beam when it hits a sample, or (false) should an
% angular distribution at the pinhole be used to model the spread of the beam.
% The model for the spread of the beam is not particuarly accurate.
allParallel = true;

% The FWHM of the density of the beam when it hits the sample (mm). For use when
% the rays are generated all parallel.
FWHM = 0.0035;

% How large is the effusive beam (proportion of the size of the main beam). Set
% to zero if the effusive beam is not to be moddeled
effuse_size = 1;

%% Pinhole plate parameters
% There are three models of the pinhole plate, varying in accuracy of the
% triangulation, 'high', 'medium', 'low', I do not believe that it should make
% much difference which is used. The manipulation of the sample will have to be
% changed if a new pinhole plate is used.
plate_accuracy = 'low';

%% Parameters for a 2d scan
% Ususally the ranges should go from -x to x. Note that these limits are in the
% coordiante system of the final image - the x axis of the final image is the
% inverse of the simulation x axis.
raster_movment2D = 0.005;
xrange = [-0.15 0.2];
zrange = [-0.2 0.2];

%% Parameters for a 1d scan
% For line scans in the y-direction be careful that the sample doesn't go
% behind the pinhole plate.
raster_movment1D = 0.1;
range1D = [-1.5 5];
Direction = 'y';

%% Sample parameters
% The sample file, include the full path
sample_fname = 'simulations/001_example/example_sample.stl';

% Sample scaling, for if the CAD model had to be made at a larger scale. 10 will
% make the model 10 times larger.
scale = 0.2;

% A string giving a brief description of the sample, for use with
% sample_type = 'custom'
sample_description = 'An example sample with some topology.';

% What type of sample to use :
%  'flat'   - A flat square (need to specify square_size)
%  'shpere' - An anlaytic sphere on a flat square surface (need to specify 
%             square_size and sphere_r)
%  'custom' - Uses the CAD model provided in the .stl file
sample_type = 'flat';

% The level of diffuse scattering for the sample, between 0 and 1, 2 gives a 
% uniform distribution. This is used for both triangulated surface and the
% analytic sphere if it is being used.
diffuse = 1;

% How close should the nearest point of the sample be to the puinhole plate, the
% defualt is 2.121 to maintain the 45o geometry. If an analytic sphere is being
% used then this is the distance between the flat surface the sphere sits on and
% the pinhole plate.
dist_to_sample = 2.121;

% The radius of the anayltic sphere (mm) (if it being included)
sphere_r = 0.05;

% If a flat sample is being used or if a sphere on a flat surface what is the 
% length of the sides of the square.
square_size = 8;


%% Output and plotting parameters

% Where to save figures/data files
% All figures and output data will be saved to this directory.
thePath = 'simulations/001_example';

% Which figures to plot
% The starting positions of the rays and the number of rays at each point
ray_starting_positions = true;
plot_density = true;

% Save simulation data to file
output_data = true;
data_fname = 'scatteringData';
saveParams = true;
paramsFile = 'scatteringParameters.txt';
% Applies only to line scans, saves in thePath with name 'data_for_plotting.csv'
save_to_text = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check some of the inputs.
% This is by no means exhaustive. Most checks are for obvious errors.
if dist_to_sample < 0.1
    error(['Cannot guarentee good results when the sample is very close to' ...
        'the pinhole plate.']);
end

if (diffuse < 0) || (diffuse > 2) || (diffuse < 2 && diffuse > 1)
    error('The given value for the scattering distribution cannot be used.');
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

if ((strcmp(sample_type, 'sphere') || strcmp(sample_type, 'flat')) && ...
        square_size <=0)
    error('The size of the flat square must be positive and non-zero.');
end

if scale <= 0
    error('The scaling of the sample model must be positive and non-zero.')
end

%% Paths to functions
addpath('plot3k', 'stlread', 'functions', 'classes', ...
        'mexFiles', 'DylanMuir-ParforProgMon-a80c9e9');

%% Path for simulation results
if ~exist(thePath, 'dir')
    mkdir(thePath)
end

addpath(thePath);

%% Generate direct beam 
% Generate the starting positions and directions of the rays
% create_starting_rays2 uses the model that I have predicted
% The source, this is the quadratic fit, ax^2, that is used as the pinhole
% angular distribution
alpha = 69.71528;

tic;

if allParallel
    % create_starting_rays creates the rays perfectly collomated with a
    % density distribution.
    [ray_pos, ray_dir, n_rays] = create_starting_rays(ray_sep, ...
        pinhole_c, pinhole_r, multipl, ray_starting_positions, thePath, ...
        FWHM, plot_density);
else
    % create_starting_rays2 creates a uniform density beam with a distribution
    % on the initial directions.
    [ray_pos, ray_dir, n_rays] = create_starting_rays2(ray_sep, pinhole_c, ...
        pinhole_r, multipl, ray_starting_positions, thePath); %#ok<*UNRCH>
end

t = toc;
fprintf('Time taken to generate rays: %f s\n', t);

%% Sample import and plotting
% Importing the sample as a TriagSurface object.
% Third argument asks if we wish to plot the surface in 3D, will save it to the
% simultion path.
switch sample_type
    case 'flat'
        sample_surface = flatSample(square_size, dist_to_sample, diffuse);
        make_sphere = 0;
        sample_description = ['A flat square sample size ' ...
            num2str(square_size) 'mm.'];
    case 'sphere'
        sample_surface = flatSample(square_size, dist_to_sample, diffuse);
        make_sphere = 1;
        sample_description = ['A single analytic sphere, radius ' ...
            num2str(sphere_r) 'mm on a flat square of ' num2str(square_size) 'mm.'];
    case 'custom'
        sample_surface = inputSample(sample_fname, diffuse, dist_to_sample, scale);
        make_sphere = 0;
end

% Do any extra manipulation of the sample here
if false
    sample_surface.moveBy([-1.27, 0, 0.05]);
    sample_surface.rotateY
end

% Plot the sample surface in 3D space, if we are using a graphical window
if feature('ShowFigureWindows')
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
        saveas(gcf, [thePath '/sample_closeUp.eps'], 'epsc');
    end
end

%% Pinhole plate import and plotting
pinhole_surface = import_plate(plate_accuracy);

% Plot if using a graphical window
if feature('ShowFigureWindows')
    sample_surface.patchPlot(true);
    pinhole_surface.patchPlot(false);
    view([-5 -5 5]);
end

%% Compile the mex files
% must include the GSL libraries
% The -ffast-math improves speed
if recompile
    mex -lgsl -lgslcblas CFLAGS="\$CFLAGS -ffast-math" mexFiles/tracingMex.c ...
        mexFiles/tracing_functions.c mexFiles/small_functions.c
    mex -lgsl -lgslcblas CFLAGS="\$CFLAGS -ffast-math" mexFiles/cosineMex.c ...
        mexFiles/small_functions.c
    mex mexFiles/binMyWayMex.c
end

%% Generate the effusive beam
% Number of rays to generate
[effuse_pos, effuse_dir, n_effusive] = makeEffuse(n_rays, effuse_size, pinhole_c);

%% Performing the simulation
switch typeScan
    case 'rectangular'
        % For a rectangular scan
        simulationData = rectangularScan(sample_surface, xrange, zrange, ...
            ray_pos, ray_dir, raster_movment2D, maxScatter, ...
            pinhole_surface, effuse_dir, effuse_pos, make_sphere, ...
            dist_to_sample, sphere_r, diffuse, thePath);
    case 'line'
        % For a line scan
        simulationData = lineScan(sample_surface, range1D, ray_pos, ray_dir, ...
            raster_movment1D, maxScatter, thePath, Direction, ...
            save_to_text, pinhole_surface, effuse_dir, effuse_pos, ...
            make_sphere, dist_to_sample, sphere_r, diffuse);
    case 'single pixel'
        % For a single pixel
        simulationData = singlePixel(sample_surface, ray_pos, ray_dir, ...
            maxScatter, thePath, save_to_text, pinhole_surface, effuse_dir, ...
            effuse_pos, make_sphere, dist_to_sample, sphere_r, diffuse);
        
        simulationData.histRays
    otherwise
        error(['Need to specify a valid type of scan: "line", ', ...
               '"rectangular", "single pixel"']);
end

%% Output data about simulation to files
if output_data
    save([thePath '/' data_fname]);
end

% Save the parameters to a text file
if saveParams 
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
    
    fprintf(fid, FORMAT1, 'Level of accuracy of pinhole plate');
    fprintf(fid, '%s\n', plate_accuracy);
    
    fprintf(fid, FORMAT1, 'Was an attempt made to model the spread of the pinhole beam?');
    if ~allParallel
        fprintf(fid, '%s\n', 'yes');
    else
        fprintf(fid, '%s\n', 'no');
        
        fprintf(fid, FORMAT1, 'The FWHM of the density of the pinhole beam');
        fprintf(fid, FORMAT2, FWHM);
    end
    
    if strcmp(typeScan, 'rectangular')
        fprintf(fid, FORMAT1, 'Number of pixels in x');
        fprintf(fid, FORMAT3, simulationData.nx_pixels);
        
        fprintf(fid, FORMAT1, 'Number of pixels in z');
        fprintf(fid, FORMAT3, simulationData.nz_pixels);
        
        fprintf(fid, FORMAT1, 'Seperation of pixels (mm)');
        fprintf(fid, FORMAT2, raster_movment2D);
        
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

