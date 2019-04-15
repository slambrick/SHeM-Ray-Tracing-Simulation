/* 
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Structues for use with the SHeM Ray Tracing Simulation.
 */

#ifndef _trace_ray_h
#define _trace_ray_h

#include "ray_tracing_structs3D.h"
#include <gsl/gsl_rng.h>

/* Function for tracing a single ray */
int32_t trace_ray_simple(Ray3D *the_ray, int *killed, int *cntr_detected, int maxScatters,
        Surface3D Sample, BackWall Plate, AnalytSphere the_sphere, gsl_rng *my_rng);

/* Trace a single ray */
int32_t trace_ray_triagPlate(Ray3D *the_ray, int *killed, int *cntr_detected, int maxScatters,
        Surface3D Sample, Surface3D Plate, AnalytSphere the_sphere,
        double backWall[], gsl_rng *my_rng);

int32_t trace_ray_simpleMulti(Ray3D *the_ray, int *killed, int *cntr_detected, 
        int maxScatters, Surface3D Sample, NBackWall Plate, AnalytSphere the_sphere,
        gsl_rng *my_rng);

#endif
