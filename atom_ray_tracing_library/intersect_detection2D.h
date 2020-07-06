/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 */

#include "ray_tracing_structs2D.h"
#include "distributions2D.h"
#include "common_helpers.h"
#include "small_functions2D.h"
#include "mtwister.h"

#ifndef _intersect_detection2D_h
#define _intersect_detection2D_h

/* Scatters all the rays stored in a Rays struct off the surface. */
void scatterRays2D(Surface2D Sample, Rays2D all_rays);

/* Traces a single ray by scattering it multiple times off a surface */
void trace_ray2D(Ray2D *the_ray, Surface2D Sample);

/* Tries to intersect a single ray with the sample surface */
int intersect2D(Ray2D *the_ray, Surface2D Sample, double intersection[2], 
                double normal[2], int *nearest_element);

#endif
