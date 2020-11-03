/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 */
#include "ray_tracing_structs2D.h"
#include "mtwister.h"

#ifndef _distributions2D_h
#define _distributions2D_h

/* Gives a ray a new direction according to the specified scattering
 * distribution. */
void new_direction2D(Ray2D *the_ray, double normal[2], int scattering, 
        double *parameters, MTRand *my_rng);

/* Creates a new direction accordin to specular reflection. */
void reflect2D(double normal[2], double init_dir[2], double new_dir[2]);

/* Creates a new direction accordin to a cosine distribution. */
void scatterCosine2D(double normal[2], double new_dir[2], MTRand *my_rng);

/* Creates a new direction accordin to a uniform distribution. */
void scatterUniform2D(double normal[2], double new_dir[2], MTRand *my_rng);

/* Creates a new direction accordin to a broad specular distribution. */
void scatterBroadSpecular2D(double normal[2], double init_dir[2], 
        double new_dir[2], double sigma_P, MTRand *my_rng);

#endif
