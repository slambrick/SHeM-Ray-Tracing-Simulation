/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 * 
 * Function declarations for scattering off different types of surface.
 */

#ifndef _scattering3D_h
#define _scattering3D_h


int scatterSphere(Ray3D *the_ray, AnalytSphere the_sphere, double *min_dist,
        double nearest_inter[3], double nearest_n[3], int *tri_hit, 
        int *which_surface);

void scatterTriag(Ray3D *the_ray, Surface3D *Sample, double *min_dist, 
        double nearest_inter[3], double nearest_n[3], int *meets, int *tri_hit,
        int *which_surface);

int backWallScatter(Ray3D *the_ray, BackWall wallPlate,  double *min_dist, 
        double nearest_inter[3], double nearest_n[3], int *meets, int *tri_hit,
        int *which_surface);

/* Scatter off a back wall with n detectors */
int multiBackWall(Ray3D *the_ray, NBackWall wallPlate, double *min_dist,
        double nearest_inter[3], double nearest_n[3], int *meets, int *tri_hit, 
        int *which_surface);

int abstractScatter(Ray3D *the_ray, AbstractHemi detector, double *min_dist, 
        double nearest_inter[3], int *meets, int *tri_hit, int *which_surface);

#endif
