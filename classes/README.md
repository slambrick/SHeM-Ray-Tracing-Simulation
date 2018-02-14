Classes for the SHeM simulation
===============================

## TriagSurface

`TriagSurface.m`

Contains a triangulated surface for use in the SHeM simulation. It is created
from the import of a binary `.stl` file, the import makes use of the code in
the `stlread` directory. Contained in the objct:

 - lists of the vertices in the triangulation
 - the indices of which vertices make up triangles (faces)
 - the unit normals to each of the faces
 - the composition of each triangle, this defines the scattering off that
   triangle, [0,1] (proportion of diffuse scattering, rest is specular) or 2
   (uniform scattering)
 - the number of triangles
 - the number of vertices.

The following methods may be of use

 - `moveBy([x,y,z])`, moves the surface by the specified amount
 - `plate_align()`, specificly for aligning the pinhole plate model into the
   correct position
 - `rotate?()`, rotates clockwise by 90 degrees about the axis ?=X,Y,Z
 - `reflectNormals`, reflects the unit normals in the surface
 - `patchPlot(new_fig, fname)`, produces a 3D plot of the surface, can specify
   if this is to be placed in a new figure (default true) or to be added to an
   already existing figure, if a file name is provided saves the figure to that
   file.

## Results classes

There are three classes that are designed to hold the results of simulation.

### SinglePixelInfo

`SinglePixelInfo.m`

Results from a single pixel simulation. Contains:

 - the number of scattering events each direct beam ray underwent
 - the final, detected, positions of each direct beam ray
 - the final, detected, directions of each direct beam ray
 - the number of rays that were forcibly stopped
 - the number of direct beam rays detected
 - the number of direct beam rays that left the simulation volume
 - the time the simulation took, in seconds
 - the number of effuse beam rays that were detected
 - the number of effuse beam rays that were forcibly stopped
 - the number of effuse beam rays that left the simulation volume

There are three methods for this class:

 - `saveDirectionsPositions(fname)`, saves the final directions and positions
   of the detected rays to a text file of the provided name
 - `single_multi_effuse()`, gives the total single, multiple and effuse detected
   rays
 - `scatteringHistogram(mmax)`, produces a histogram of the number of detected
   rays that underwent each number of sample scattering events, returns the
   histogram as an array and produces a plot.

### LineScanInfo

`LineInfo.m`

Results from a line scan simulation. Contains:

 - the direction of the line scan, 'x','y','z'
 - a list of the sample positionsin the direction of the scan
 - the range of the scan [min, max]
 - the number of pixels in the scan
 - an array of histograms of the number of rays that underwent 1,2,3,4...
   sample scattering events for each pixel, for the direct beam
 - the number of detected rays from single scattering for each pixel, from the direct beam
 - the number of detected rays from multiple scattering for each pixel, from the direct beam
 - the number of rays that had to be forcibly stopped
 - the number of detected rays from single scattering for each pixel, from the effuse beam
 - the number of rays from the effuse beam that had to be forcibly stopped
 - the distance the sample was moved between pixels
 - the number of rays for each pixel
 - the time the simulation took, in seconds
 - an estimate of the time the simulation would take, in seconds.

There are three methods for this class:

 - `saveText(fname)`, saves the final positions and directions of the direct
   beam rays to a text file of the specified name
 - `totalLine()`, gives the total (sum of all contributions) line scan
 - `producePlots()`, produces two plots of the line scan, number of detected
   rays as a function of scan position. One seperates the different
   contributions, the other is the total.

### RectangularInfo

`RectangleInfo.m`

Results from a 2D simulation. Contains:

 - an array of histograms of the number of rays that underwent 1,2,3,4...
   sample scatterng events for each pixel, for the direct beam
 - a matrix of the total number of rays detected in each pixel for the direct
   beam
 - a matric of the total number of rays detected in each pixel for the effuse
   beam
 - the number of pixels in the *x* direction, the 'horizontal' direction
 - the number of pixels in the *z* direction, ther 'vertical' direction
 - the total number of pixels
 - the number of direct beam rays in each pixel
 - a TriagSurface object of the sample
 - the range of the scan [min, max] in the x direction
 - the range of the scan in the z direction
 - the movement of the sample between pixels
 - the time the simulation took, in seconds
 - an estimate of the time the simulation would take, in seconds

Ther are 9 public methods and 1 private method for this class. All of them
create a 2d plot of the results in some fashion. One method
`produceImages(path)` produces a series of images and saves them to the provided
directory. Images of the multiple scattering, single scattering, effuse, total
direct and total contributions are saved as `.png` images, with axes. It also
produces and saves to a `.eps` file a contour plot of the image.

The rest of the methods produce a single visualisation of the data. They have
the generic form

```matlab
imageX(scale, specifyScale, limX, limY)```

X is the type of image to produce. `scale` is how to produce the grayscale:
 - 'zeroed', makes black be 0 counts, white is the pixel with the most counts
 - 'auto', makes the pixel with the fewest counts black, white is the pixel
   with the most counts
 - 'removeZero', makes black the smallest number of detected rays that isn't
   zero
 - 'manual', manually provide how many counts correspond to black and white

If no `scale` is provided then 'auto' is chosen. Note that `imageBigger` can
only use 'auto' and 'zeroed'. If 'manual' is specified then `specifyScale` must
also be provided, this is a two element vector [min, max]. The x and y limts of
the images may be specified, this requires a `specifyScale` to be specifed even
if it isn't used. `limX` and `limY` give the limits of the horizontal and
vertical axes respectivly.

`imageExtraEffuse(factor, scale,...)` allows the effuse contribution to be
increased or decreased by a provided factor (multiplies that contribution by the
factor.

`contourImage(n, limX, limY)` produces a contour plot of the total direct beam
contribution. `n` is the number of contours to plot, it defaults to 25. `n`
must be specified if `limX` and `limY` are to be specified.
