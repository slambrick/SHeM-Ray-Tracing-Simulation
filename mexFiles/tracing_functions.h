/* 
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Function declerations for the larger ray tracing function for the SHeM
 * Simulation.
 */

#ifndef _tracing_functions_h
#define _tracing_functions_h

//#include <gsl/gsl_rng.h>
#include "ray_tracing_structs3D.h"

/* 
 * Finds the intersection, normal at the point of intersection and distance to 
 * the intersection between the ray and an analytic sphere. returns 0 if the 
 * ray does not intersect the sphere. 
 */
int scatterOffSurface(Ray3D *the_ray, Surface3D *Sample, AnalytSphere 
    the_sphere);

/* 
 * Finds the intersection, normal at the point of intersection and distance to 
 * the intersection between the ray and a triangulated surface. 
 */
int scatterPinholeSurface(Ray3D *the_ray, Surface3D *Plate, double backWall[]);

/* 
 * Scatters a ray off two triangulared surfaces, and an analytic sphere if 
 * desired. 
 */
int scatterSurfaces(Ray3D *the_ray, Surface3D *Sample, Surface3D *Plate, 
        AnalytSphere the_sphere, double BackWall[]);

/* 
 * Scatters a ray off a triangulated surface, and a simple model of the pinhole plate
 * as a surface.
 */
int scatterSimpleSurfaces(Ray3D *the_ray, Surface3D *Sample, BackWall Plate, 
        AnalytSphere the_sphere);

/* 
 * Scatters a ray off a triangulated surface, and a simple model of the pinhole plate
 * with multiple detector apertures.
 */
int scatterSimpleMulti(Ray3D *the_ray, Surface3D *Sample, NBackWall Plate, 
        AnalytSphere the_sphere, int *detector);

/* 
 * Scatters the ray off a sample and a attempts dtection on a hemisphere with abstract
 * detector apertures placed on it.
 */
int scatterAbstractSurfaces(Ray3D *the_ray, Surface3D *Sample, AbstractHemi Plate,
        AnalytSphere the_sphere) ;

#endif
