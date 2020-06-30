/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 * 
 * Function declarations for the small functions used in the SHeM ray tracing 
 * simulation.
 */

#ifndef _scattering_functions_h
#define _scattering_functions_h

#include "mtwister.h"
#include "ray_tracing_structs3D.h"

void new_direction3D(Ray3D *the_ray, double normal[3], double scattering, 
    double parameters, MTRand *myrng);

void broadSpecular3D(double normal[3], double init_dir[3], double new_dir[3],
        double sigma, MTRand *myrng);

void cosineScatter3D(double normal[3], double new_dir[3], MTRand *myrng);

void cosineSpecularScatter3D(double normal[3], double initial_dir[3], 
        double new_dir[3], MTRand *myrng) ;

void uniformScatter3D(double normal[3], double new_dir[3], MTRand *myrng);

void reflect3D(double normal[3], double init_dir[3], double new_dir[3]);

#endif
