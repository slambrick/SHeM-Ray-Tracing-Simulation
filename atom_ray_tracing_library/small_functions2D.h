/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 */
#include "ray_tracing_structs2D.h"

#ifndef _small_functions2D_h
#define _small_functions2D_h

/* Dot product between two 2 element vectors. */
double dot2(double a[2], double b[2]);

/* The square of the norm of a 2 element vector. */
void norm2(double vect[2], double* result);

/* Normalises a 2 element vector. */
void normalise2(double vect[2]);

/* Prints out the vertices and normals of the sample surface. */
void print_surface2D(Surface2D sample);

/* Gets the nth element of a 2D Surface struct and put the information in the
 * provided arrays. */
void get_element2D(Surface2D Sample, int n, double v1[2], double v2[2], 
                   double nn[2]);

/* Puts the rays defined by the arrays of positions and directions into an array
 * array of Ray structs. */
void compse_rays2D(Ray2D *rays, double ray_pos[], double ray_dir[], int nrays);

/* Prints a ray position and direction */
void print_ray2D(Ray2D ray);

/* Prints all the 2D rays in the Rays struct. */
void print_all_rays2D(Rays2D all_rays);

/* Puts the directions of the rays into a single matrix for outputting. */
void get_directions2D(Rays2D all_rays, double final_dir[]);

/* Puts the positions of the rays into a single matrix for outputting. */
void get_positions2D(Rays2D all_rays, double final_pos[]);

#endif
