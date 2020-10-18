/*
 * Copyright (c) 2018-20, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Function declarations for the small functions used in the SHeM ray tracing
 * simulation.
 */

#ifndef _small_functions3D_h
#define _small_functions3D_h

/* Adds a 3 element array to another 3 element array multiplied by a scalar
 * (propagates an array). */
void propagate(const double init[3], const double direc[3], double a, double result[3]);

/* Calculates the dot product of two 3 element double vectors. */
void dot(const double a[3], const double b[3], double* result);

/* Calculates the cross product of two 3-vectors and writes it to c */
void cross(const double a[3], const double b[3], double c[3]);

/* Returns the square of the norm of a three vector. */
void norm2(const double vect[3], double* result);

/* Normalises a three vector */
void normalise(double vect[3]);

/* reflect direction through normal */
void reflect3D(const double normal[3], const double init_dir[3], double new_dir[3]);

/* find two (unit) directions perpendicular to a given unit vector,
 * and write them to v1 and v2. */
void perpendicular_plane(const double n[3], double v1[3], double v2[3]);

/* Solves a 3D matrix equation Au=v. */
void solve3x3(double A[3][3], double u[], double v[], double epsilon, int* success);

#endif
