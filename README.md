SHeM ray tracing simulation
===========================

## Principles

The code here is used to perform a simulation of the Cambridge scanning helium
microscope (SHeM), though it is easily adapted to other similar atomic
microscopes, and with more difficulty may be applied to other forms of
microscopy. The basic idea consists of modelling the atoms in the SHeM as rays
that travel in straight lines and scatter randomly according to a predefined
distribution when they hit a surface. The surfaces are defined by a
traiangulated mesh. The brightness of a pixel is determined by the number of
rays that reach a surface that is defined to be the detector.

To gain a fuller understanding of the principles of the simulation and its
capabilities it is recommended to read the the associated publication (coming
soon... ). The design of the Cambridge SHeM is detailed in
[A design for a pinhole scanning helium microscope](https://doi.org/10.1016/j.nimb.2014.06.028)
M. Barr *et al.* 2014.

---

## The simulation

### Structure of the simulation

A main `.m` script `performScan.m` contains the specification of the parameters
for the simulation and when run calls functions that perform the simulation.
Functions exist for each of the three types of simulation: a *single pixel*
simulation, where only one pixel is taken with the sample in a fixed position;
a *line scan*, where a line of pixels are taken with the sample moving in a
line; or a *rectangular scan*, which takes pixels in a rectangular region
moving the sample appropriatly. Each of these functions then calls a wrapper
function that itself calls the `.c` MEX files that perform the core ray tracing.

A series of other functions and MATLAB classes assist in the setup of the
simulation and in the processin of the output. Custom classes contain the
results of simulations, one for each type of scan. More information is retained
in the case of a single pixel simulation than in the other two options. The only
parameter that is not totally specified within `performScan.m` is the 3D model
of the sample. These are placed within the simulation results directory and the
file name is then specified in `performScan.m`.

### Subdirectories

Each subdirectory in the main simulation directory (the current one) contains
components of the simulation. All scripts to be run are found in the main
directory. `.m` code is found in the *functions* and *classes* folders, the
files therin contain what would be expected. `README.md` files in those folders
detail what those functions and classes do. The `.c` code that makes up the MEX
files that perform the core ray tracing are found in the *mexFiles* folder.

The folder *simulations* is where simulation results are placed. Each simulation
performed should go in its own folder within simulations, e.g. *001_example*.
The folder *pinholePlates* contains binary `.stl` files containing models of the
Cambridge SHeM pinhole plate, the orginal design along with three accuracies
of a simplified model to be used in the simulation. The simplest `.stl` file
should suffice for the simulations.

The following contain code that does not belong to
the main author:

 - *DylanMuir-ParforProgMon-a80c9e9*, contains code for producing a graphical
   progress bar that works with parfor loop, makes use of java.
 - *stlread*, contains code for importing binary `.stl` files.
 - *plot3k*, contains code for producing 3D scatter plots.

Each of these folders contains the licence for the respective codes.

---

## Software prerequisites for running the simulation

### MATLAB

The simulation is written largely in MATLAB `.m` code, it has been successfully
run on the following versions of MATLAB

 - r2014b
 - r2015b
 - r2017b.

The simulation makes use of the parallel computing toolbox and the image
processing toolbox, however, dependence on these can easily be removed and the
simulation will still work (though slower and there will be no method for
displaying the results for 2D simulations).

### Compiling MEX files

The core ray tracing element of the simulation is performed using C code that
must be compiled for the particular system being used. The compilation is achieved
using MEX. The GNU scientific library (GSL) is also used, in particular the
random number generator library, `gsl/gsl_rng.h`. In order to run the simulation
the GSL for C must be installed, for compilation the development library for the
GSL must also be installed (the development packes are needed for compiling code
that uses the GSL). In Ubuntu 17.10 the development packages is called `libgsl-dev`
and the non-development packages are dependent on it.

Compilation has been tested successfully using:

 - gcc 7.2.0 on Ubuntu 17.10
 - gcc 4.8.4 on Ubuntu 14.04 LTS
 - gcc 4.7.2 on Scientific Linux 7.3

It should be noted that only only some gcc versions are officially supported with
MATLAB and MEX. Most of the versions tested are not officially supported.

### Hardware requirements

The simulation is CPU intensive, though will run on any CPU that can run MATLAB.
RAM usage will vary significantly depending on the simulation being run.
Increasing the complexity of the sample or increasing the number of rays used
per pixel will increase the RAM usage. As a ballpark I recommend there being 2GB
of RAM per core being used to run the simulation.

### GNU Octave

The simulation cannot be run in GNU Octave 'as is', there are a number of features
used that are not available.

 - The use of `parfor` is not compatible (there are other
ways of parallelising code in Octave).
 - The Octave Forge 'statistics' package is
needed.
 - The simulation only tries to display figures when run from a graphical
window, the function `feature()` that is used to determine whether a graphical
window is being used does not exist in Octave, those sections of code must be
edited.
 - The package `liboctave-dev` needs to be install (Debian based systems)
for compilation of the MEX files.
 - The option `CFLAGS="\$CFLAGS -ffast-maths` does not work when used in Octave.
 - The code that generates the progress bar for the simulation need to be disabled.
 - A number of warnings are shown when the MEX files are executed.
 - The results classes are set to be immutable' using the code `(SetAccess = immutable)`, this needs to be removed.
 - AN alternative method of saving the results must be used `save('results.mat')` doesn't work.
 - The simulation runs slower than it would using MATLAB.

---

## Simulation files

### Scripts

In the main directory in addition to the main simulation script there are two
other scripts, these are used for testing scattering distribution and how the
pinhole beam spreads when a directional distribution is placed on the initial
rays.

`performScan.m` is the main simulation script, it is run from within MATLAB
to perform the SHeM simulation.

`test_cosineDist.m` tests the sampling of the cosine scattering distribution
as it is done from within the simulation. The script is self explanitory and
uses its own MEX gateway function. It can be adapted to test the sampling
of a different scattering function by modification of the gateway function.

`beam_spread.m` explores the spread of the pinhole beam as it travels away from
the pinhole, given a directional spread from the pinhole. In the SHeM
simulation thera are two basic models of the pinhole beam. The one that is in
use in the `performScan.m` as it is proivided uses a perfectly collimated beam
that has been given a density, a 2D Gaussian is used. The other model is to
assume a unifrom density across the pinhole and to vary the direction of rays
that come out of the pinhole. Details of the moddeling of the pinhole beam in
such a way are given below.

### Classes

Four classes are used in the SHeM simulation. Three of them contain the results
of simulations, the fourth, `TriagSurface`, contains the triangulation of a
surface. This class is created by the inport function `inputSample.m` and is
used for both the sample and the pinhole plate. It is possible to plot such an
object of `TriagSurface` with `exampleSurface.patchPlot()`.

The following classes are used to contain simulation results:

 - `SinglePixelInfo`
 - `LineScanInfo`
 - `RectangleInfo`

### `.m` functions

A number of functions written in `.m` code are to be found in *functions*, these
functions are used to set up the simulation, then provided useful wrappers of
MEX functions to perform the main types of simulation easily, and finally do
analysis of the results and save the data.

### MEX files

There are a series of MEX files that perform parts of the simulation, and some
small `.c` libraries. The core ray tracing is performed using the function/file
`tracingMex.c`, two other MEX functions, `binMyWayMex.c` and `cosineMex.c`
perform smaller parts of the simulation. More details of these functions can be
found in the `README.md` in *MexFiles*.

---

## Running the simulation

### Setting up the parameters

All the parameters for the simulation are specified in the script
`performScan.m`. The script is separated into three main sections, the first
should be edited to specify all the parameters, the second calls a series of
functions to perform the simulation, and the third outputs the results to a
MATLAB `.mat` file and some of the parameters. The third section can easily be
edited if different parameters or data are wished to be saved, editing the
second part should be approached more carefully. For multipixel simulations the
program makes use of the parallel computing toolbox, it runs the simulations for
separate pixels on separate cores, if this functionality is not desired simply
change the `parfor` loops in the function `rectangularScan.m` and `lineScan.m`
to `for` loops and remove the lines that create and close the MATLAB parallel
pools.

A separate file, `how_to_run.md`, lists all the parameters and how to specify
them properly.

It should be noted that when running any multipixel scan a time estimate is
printed out and the user is asked if they wish to proceed, this is because
simulations can take many hours, or days, if a complicated simulation is being
done. The time estimate is very much a rough 'order of magnitude' estimation.

### Interpreteting and processing results

After a simulation is run an object called `simulationData` is created and
saved to a `.mat` file in the specified simulation directory. This object will be
of one of the three classes `SinglePixelInfo`, `LineInfo`, or `RectangleInfo`
depending on if a single pixel simulation, line scan, or a full 2D scan was
undertaken. This object will have all the results in a reasonably logical
fashion. In the cases of a single pixel or line scan simulation some of the data
will also be saved to a text `.csv` file in the same directory, if analysis or
plotting wants to be performed in external software.

For the `SinglePixelInfo` and `LineScanInfo` classes there methods that can
produce some basic plots of the results, neither of them automatically save the
plots. For the line scan `producePlots` creates two graphs of the number of
detected rays as a function of the line scan position, one with just the total
number of rays and another seperating the results into single scattering,
multiple scattering, and effuse beam contributions. For a single pixel
simulation `scatteringHistogram` produces a histogram of the number of detected
rays that had undergone 1,2,3... number of sample scattering events.

For 2D simulations a series of images are automatically produced and saved in
the simulation directory. One of the single scattering contribution, one of the
multiple scattering contribution, one of the effuse contribution, one of the
total without the effuse beam, and one of the total; a contour plot of the image
is also produced. All of the images are saved as `.png` files with axes in mm,
and the contour plot is saved as an `.eps` file. There are a series of methods
for `RectangleInfo` that allow different types of image to be easily
constructed.

There is an `.R` script in the simulations directory that produces plots for line
scans, running this script specifying where the results files are and which
which direction the line scan was in produces a series of plots in the directory
of simulation results. The script makes use of packages in the 'tidyverse'. MATLAB also produces line plots.

---

## Spreading of the pinhole beam

The script `beam_spread.m` samples the angular distribution of atoms predicted
to be at the pinhole assuming that the skimmer in the SHeM is the source. It
then propagates rays given those directions a distance that may be specified
and plots their final position overlayed with the pinhole.

---

## Licence notice

Copyright (c) 2018, Sam Lambrick. 
All rights reserved.
The SHeM Ray Tracing Simulation is subject to the GNU General Public Licence 
version 3.0-or-later. The licence is found in `LICENCE` in the main directory
of this repository or at <https://www.gnu.org/licenses/gpl.html>.

---

## Credit

It is requested that if use of this work results in an academic publication that
it is appropriately cited and that if the ideas that led to the creation of the
simulation are used then the associated publication should be cited (coming
soon... ).
