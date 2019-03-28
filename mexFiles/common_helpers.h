/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 */

#ifndef _common_helpers_h
#define _common_helpers_h

#include <gsl/gsl_rng.h>

/* Sets up the GNU/SL random number generator, must be called before any 
 * function that uses random numbers. */
gsl_rng* setupGSL(void);

/* Prints a 2 or three element vector */
void print1D_double(double *vect, int dim);

/* Prints a 2 or three element vector */
void print1D_int(int *vect, int dim);

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double matrix[3][3]);

#endif
