Functions for the SHeM simulation
=================================

The `functions` directory contains `.m` functions that are used as part of the
SHeM Ray Tracing Simulation. This file lists the functions in this directory
and goes through what they are used for.

#### `binMyWay`

A function for producing a specific type of histogram, it is used to sort the
detected rays by the number of scattering events they have undergone. Wrapper
for a MEX function.

#### `combineRectangle`

Combines two sets of simulation data for two equivalent 2D scans. Equivalent
means that they are of the same sample with the same parameters, the number of
rays is an exception. This can be used when a particularly computationally
expensive simulation is needed, the simulation can be split up and performed
on different computers, each of the same area of sample but with too few rays.

#### `composition_sample`

Creates a sample that shows how scattering can be made to be different from
different parts of a sample.

#### `create_starting_rays`

Creates the starting rays for a parallel beam with a Gaussian profile.

#### `create_starting_rasy2`

Creates the starting rays for a beam that uses a directional distribution.

#### `flatSample`

Creates a flat sample of the specified size. Used to make the surface the
analytic sphere sits on.

#### `import_plate`

Specifically for importing the pinhole plate.

#### `inputSample`

Imports a sample from a `.stl` file and makes a `TriagSurface` object out of it.

#### `lineScan`

Performs a 1D scan across the sample in one of three directions.

#### `makeEffuse`

Creates the effuse beam.

#### `random_dir`

Generates the directions for the rays for use with `create_starting_rays2`.

#### `rectangularScan`

Performs a rectangular 2D scan to produce a 2D simulated image of the sample.

#### `singlePixel`

Performs a simulation with the sample in a fixed position.

#### `time_estimate`

Makes an estimate of the time the simulation will take using the number of rays
and the complexity of the sample. Prints out to the screen and asks the user
if they want to continue.

#### `trace_rays`

Wrapper function for the C ray tracing functions. This should be called by one
of the three functions `lineScan`, `rectangularScan`, or `SinglePixel` and not
directly.