SHeM ray tracing simulation
===========================

This repository contains the SHeM ray tracing simulation that can generate simulated
images (and similar data) for scanning helium microscopy. The basic principle 
of the simulation is that the paths of atoms follow straight lines, only changing 
direction upon a scattering event with a surface. The simulation can include more or
less realistic models of the local sample environment and can be used either with
simple test samples or with your own complex samples drawn in CAD.

## Getting/running the code

See the wiki for detailed information on the simulation and instructions on how 
to run it, only some brief notes are given below. 

1. Download the code or clone the repository
2. Ensure you have the prerequisit software (MATLAB + supported complier, see the wiki for details)
3. Complie the MEX code `mexCompile` (this will be done automatically when you run your first simulation
4. Set the parameters for your simulation in the `ray_tracing_parametes.txt` file
5. Run the main simulation script, `performScan.m`

## Licence notice

Copyright (c) 2018, Sam Lambrick. 
All rights reserved.
The SHeM Ray Tracing Simulation is subject to the GNU General Public Licence 
version 3.0-or-later. The licence is found in `LICENCE` in the main directory
of this repository or at <https://www.gnu.org/licenses/gpl.html>.

---

## Credit

It is requested that if use of this work results in an academic publication that
it is appropriately cited (DOI: [10.5281/zenodo.1228079](https://doi.org/10.5281/zenodo.1228079)) and that if the ideas that led to the creation of the
simulation are used then the associated publication should be cited: [A ray tracing method for predicting contrast in neutral atom beam imaging](https://doi.org/10.1016/j.micron.2018.06.014).
