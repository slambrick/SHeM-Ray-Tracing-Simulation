Detail of the parameters and running SHeM simulations
=====================================================

This file walks through the parameters of the simulation.

## General parameters

Here some general simulation parameters are specified. The maximum number of sample
scatters per ray can be specified. If the simulation is being run on a computer/
operating system for the first time then it is important to set `recompile` to
be `true`. If the MEX files are not compiled for the system being used then the simulation can't run. The type of scan being performed is also specified.

```Matlab
% The maximum number of sample scatters per ray. There is a hard-coded total
% maximum number of scattering events of 1000 (sample and pinhole plate). Making
% this unnecessarily large will increase the memory requirements of the
% simulation.
maxScatter = 20;

% Type of scan 'line', 'rectangular', or 'single pixel'
typeScan = 'rectangular';

% Recompile mex files?
% Required if using on a new computer or if changes to .c files have been made.
recompile = true;
```

---

## Beam/source parameters

The nature of the pinhole beam must be specified. First define the centre of the
pinhole in the coordinate system of the simulation, in general keep the y and z
coordinate zero. The the total extent of the pinhole must be specified, there are
two models for the pinhole, if the model that uses the skimmer as the source is
being used then the dimension specified here is the true radius of the pinhole,
if a Gaussian beam profile (default) is being used it is the radius of the region
over which the beam profile is being modelled (larger than the FWHM of the beam).

The number of rays used is specified by a combination of: the separation of the
rays in the pinhole, they are created in a square grid with the specified
separation; the extent of the beam, defined by `pinhole_r`; the maximum number
of rays at each location, the total number at each location if the skimmer-source
model is being used; the FWHM of the Gaussian beam if a beam profile is used.

Finally the relative size of the effuse beam must be specified, this is relative
to the size of the direct beam, i.e. `effuse_size=1` means have the same number of
rays in both. The effuse beam follows a cosine distribution from the pinhole
centred on the normal to the pinhole plate.

```Matlab

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
```

---

## Pinhole plate parameters

The accuracy of the pinhole plate must be specified, 'low' should almost always
be used.

```Matlab
% There are three models of the pinhole plate, varying in accuracy of the
% triangulation, 'high', 'medium', 'low', I do not believe that it should make
% much difference which is used. The manipulation of the sample will have to be
% changed if a new pinhole plate is used.
plate_accuracy = 'low';
```
---

## Parameters for a 1D and 2s scans

The movement between pixels and the range of the scans must be specifed. For a
line scan the direction must also be specified. 'x' is the horizontal axis, 'z'
is the vertical axis (which is the y axis when an image is produced), and 'y' is
corresponds to diagonally moving the sample away from the pinhole plate, with
the same region of sample in the beam.

```Matlab
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
```

---

## Sample parameters

There are three types of sample that can be used by default. A flat square,
a analytically modelled sphere sat on a flat square, or a custom sample
defined in a binary `.stl` file. Which type is specified by `sample_type`.

In the case of a custom sample then the name of the `.stl` file including the
full path to the file must be specified as must any scaling of the sample and
a short description of the sample. A scaling of 0.1 makes the imported sample
10 times bigger and 10 makes the imported sample 10 times smaller.

In the case of a flat sample the size of the square needs to be specified.

In the case of a sphere the radius of the sphere and the size of the square it
sits on must be specified.

In all cases the scattering off of the sample must be specified:
 - 2 = uniform scattering
 - 1 = normally centred cosine
 - 0 = completely specular
 - between 0 and 1 = a combination of specular and cosine (0.8 means 20% specular)

The distance from the pinhole plate to the sample must also be specified. For the
cases of a custom sample or a flat sample it is the minimum distance for any point
on the sample. For the case of an analytic sphere it is the distance between the
pinhole plate and the flat square the sphere sits on.

```Matlab
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
%  'sphere' - An anlaytic sphere on a flat square surface (need to specify 
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
```

---

## Output and plotting parameters

Where to save data, what to name the files that contain the data and which
figures to plot. `output_data` specified if to save the data to a `.mat` file
named by `data_fname`, `saveParams` specifies if the parameters should be saved
to a text file named by `paramsFile`, `save_to_text` specified if data from
line scans simulations should be saved to a text file (for plotting 
or analysis in other software. The starting positions of the rays in the pinhole 
can be plotted, as can the density of the rays if the Gaussian profile is being 
used.

```Matlab
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
```
