Mex files in SHeM simulation
============================

The core of the SHeM simulation is performed in C. All the C/MEX files used are
in this folder (mexFiles). Fiddlling with these function without proper
knowledge can may break the simulation, just a warning, the simulation can be
adapted to many uses without having to go down here.

Introduction of new scattering distributions requires the addition of functions
and modification of excisting code. It is recommended that any new scattering
functions are added to `small_functions.c`.

## Prerequisits

The files in this folder are written to be complied via mex inside MATLAB. They
make use of the `mex.h`, `math.h`, `stdio.h`, and `time.h` libraries as well as
the standard library  `stdlib.h`. The GNU scientific library (GSL) is also used,
in particular the random number generator library, `gsl/gsl_rng.h`. In order to
run the simulation the GSL for C must be installed, for compilation the
development library for the GSL must also be installed. In Ubuntu 17.10 the
development packages is called `libgsl-dev`.

Compilation has been tested successfully using:

 - gcc 7.2.0 on Ubuntu 17.10
 - gcc 4.8.4 on Ubuntu 14.04 LTS
 - gcc 4.7.2 on Scientific Linux 7.3

It should be noted that only only some gcc versions are officially supported with
MATLAB and MEX. Most of the versions tested are not officially supported.

## The C mex files used in the imaging simulation.

`tracingMex.c` is the mex file that is to be called by the MATLAB wrapper
function and performs the core ray tracing computation for a single pixel. It is
a complicated function that takes many inputs and gives many outputs. More data
is outputted than is used in any of the multiple pixel simulations, look at the
ouptuts of the wrapper function to see what the simulation can tell you.

`cosineMex.c` is the mex file for sampling a cosine distribution, it is used to
test the function that samples the distribution and also to generated the effuse
beam in the SHeM simulation.

`binMyWayMex.c` is a small function that bins data to produce a histogram, I
wrote it as I was having trouble getting the MATLAB functions to do exactly what
I wanted. It bins an integer array (of at least length three) into a number of
bins specified, with the first bin being 1 and going up in integer values. It
should be used only for its specific purpose, use histogram (or hist... ?)
instead for your binning.

## C libraries for the imaging simulation

`tracing_functions.c` contains larger functions that perform the tracing of a
single ray off of the surfaces in the simulation. Two main functions are used
to scatter the ray off of either all the surfaces or the surfaces without the
pinhole plate. Other functions calculate the intersection between surfaces and
rays.

`small_functions.c` contains a number of smaller functions that are used in the
main algorithms in `tracing_functions.c` and `tracingMex.c`. There are also
functions that are useful for debugging, such as printing functions, and
functions no longer used but are generally useful for the ray tracing as it is
performed here.
