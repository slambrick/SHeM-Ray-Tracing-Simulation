/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * Functions that trace a single ray through multiple collisions
 * until it is either detected or dead.
 */

#ifndef _trace_ray_h
#define _trace_ray_h

#include "ray_tracing_core3D.h"
#include "mtwister.h"
#include <stdint-gcc.h>


void trace_ray_simple_multi(Ray3D *the_ray, int maxScatters, Surface3D sample,
        NBackWall plate, AnalytSphere the_sphere, MTRand * const myrng);

void trace_ray_triag_plate(Ray3D * the_ray, int maxScatters, Surface3D sample,
        Surface3D plate, AnalytSphere the_sphere, double const backWall[], MTRand * const myrng);

void trace_ray_just_sample(Ray3D * the_ray, int * const killed, int maxScatters,
        Surface3D sample, AnalytSphere the_sphere, MTRand * const myrng);

#endif
