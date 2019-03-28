/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 * 
 * Function declarations for the small functions used in the SHeM ray tracing 
 * simulation.
 */

#ifndef _small_functions_h
#define _small_functions_h

/* Adds a 3 element array to another 3 element array multiplied by a scalar 
 * (propagates an array). */
void propogate(double init[3], double direc[3], double a, double result[3]);

/* Calculates the dot product of two 3 element double vectors. */
double dot(double a[3], double b[3]);

/* Returns the square of the norm of a three vector. */
double norm2(double vect[3]);

/* Normalises a three vector */
void normalise(double vect[3]);

/* Solves a 3D matrix equation Au=v. */
int solve3x3(double A[3][3], double u[], double v[], double epsilon);

#endif
