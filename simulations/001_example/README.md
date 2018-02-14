Example sample for the simulation
=================================

In this directory is an example sample for the simulation, `example_sample.stl`
and the Autodesk Inventor part file that it was created from. To import the
sample a scaling of 0.2 should be used. Then a scan area from -0.2 to 0.2mm in
z and -0.15 to 0.2 in x can be used. Make `raster_movment2D` relativly large,
~0.005 or the simulation will take a very long time to run. The example figures
do not included the effuse beam.

To set the number of rays the following set of parameters works:
 - pinhole radius = 0.005mm
 - max number of rays at each point = 40
 - ray separation in the pinhole = 0.0001mm
 - FWHM of the pinhole beam (Gaussian) = 0.0035mm
