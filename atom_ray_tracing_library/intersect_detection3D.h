/*
 * Copyright (c) 2018-20, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Function declarations for scattering off different types of surface.
 */

#ifndef _intersect_detection3D_h
#define _intersect_detection3D_h

#include <stdbool.h>

#include "ray_tracing_core3D.h"

void scatterSphere(Ray3D * the_ray, AnalytSphere const * const the_sphere, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], int * const tri_hit,
        int * const which_surface, int * const meets_sphere);

void scatterTriag(Ray3D * the_ray, Surface3D const * const Sample, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], int * meets, int * const tri_hit,
        int * const which_surface);

void multiBackWall(Ray3D * the_ray, NBackWall const * const wallPlate, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], int * const meets, int * const tri_hit,
        int * const which_surface, int * const which_detector);

//void abstractScatter(Ray3D * the_ray, AbstractHemi const * detector, double * min_dist,
//        double nearest_inter[3], bool * meets, int * tri_hit, int * which_surface,
//		int *which_detector);

#endif
