Mex files in SHeM simulation
============================

The ray tracing element of the simulation in performed in C and mex is used to 
pass information from Matlab/Octave to C. The C code is split into three parts,
first the mex functions that can be called from Matlab/Octave, then files
that genereate the particular simulations of a ray in some geometry, and finally
a library for atom ray tracing.

## Prerequisits

Compilation has been tested successfully using:

 - Ubuntu 18.04  
 - gcc 7.2.0 on Ubuntu 17.10
 - gcc 4.8.4 on Ubuntu 14.04 LTS
 - gcc 4.7.2 on Scientific Linux 7.3

It should be noted that only only some gcc versions are officially supported with
MATLAB and MEX. Most of the versions tested are not officially supported.

## The mex functions

## Simulation functions

## The C atom ray tracing library
