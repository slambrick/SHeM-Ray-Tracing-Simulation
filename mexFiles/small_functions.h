/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 * 
 * Function declarations for the small functions used in the SHeM ray tracing 
 * simulation.
 */

#ifndef SMALL_FUNCTIONS_H
#define SMALL_FUNCTIONS_H


/* Prints out a 3 element 1D double array passed as an argument. */
void print1D(double vect[3]);

/* Prints out a 3 element 1D int array passed as an argument. */
void print1Dint(int vect[3]);

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double matrix[3][3]);

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

/* Specuarly reflects a ray, given the initial direction and the surface 
 * normal */
void reflect(double normal[3], double init_dir[3], double new_dir[3]);

/* Seeds the C standard library random number generator. */
void seederForScattering(void);

/* Seeds the GNU/SL random number generator */
gsl_rng* setupGSL(void);

/* Generates a random normalized direction according to a cosine distribution
 * about the provided normal and stores the result in the provided array. */
void cosineScatter(double normal[], double new_dir[], gsl_rng *myrng);

/* Generates a random normalized direction according to a uniform distribution
 * about the provided normal and stores the result in the provided array. */
void uniformScatter(double normal[], double new_dir[], gsl_rng *myrng);

#endif
