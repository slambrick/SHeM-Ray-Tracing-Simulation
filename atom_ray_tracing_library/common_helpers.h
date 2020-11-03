/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */
 
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif


#ifndef _common_helpers_h
#define _common_helpers_h

#include "mtwister.h"


/*
 * Linearise [row][column] coordinates in an array of coordinates
 * such as a 3-column n-row matrix of vertices of a Surface3D
 */
void lin(int row, int col, int * const ind);

/* Prints a 2 or three element vector */
void print1D_double(double const * const vect, int dim);

/* Prints a 2 or three element vector */
void print1D_int(int const * const vect, int dim);

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double const matrix[3][3]);

void gaussian_random(double mu, double sigma, double Z[2], MTRand * const myrng);

void gaussian_random_tail(double mu, double sigma, double cutoff, MTRand * const myrand,
		double * const rand1);

void gen_random_int(int max, MTRand * const myrand, int * const randint);

#endif
