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

void scatterSphere(Ray3D * the_ray, AnalytSphere * the_sphere, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const tri_hit,
        int * const which_surface, bool * const meets_sphere);

void scatterPlane(Ray3D * the_ray, Plane plane, double * const min_dist,
		double nearest_inter[3], double neares_n[3], double nearest_b[6], int * const meets,
		int * const which_surface);

void scatterCircle(Ray3D * the_ray, Circle * the_circle, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const tri_hit,
        int * const which_surface, bool * const meets_circle);

void scatterTriag(Ray3D * the_ray, Surface3D * Sample, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * meets, int * const tri_hit,
        int * const which_surface);

void scatterSample(Ray3D * the_ray, Sample overall_sample, double * min_dist, double * nearest_inter, 
                   double * nearest_n, double * nearest_b, int * meets, int * tri_hit, int * which_surface);

void multiBackWall(Ray3D * the_ray, NBackWall wallPlate, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const meets, int * const tri_hit,
        int * const which_surface, int * const which_aperture);

void abstractScatter(Ray3D * the_ray, AbstractHemi const * plate, double * const min_dist, 
        int * const tri_hit, int * const which_surface, int * const which_aperture);

#endif
