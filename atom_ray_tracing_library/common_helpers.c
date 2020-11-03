/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Contains small helper functions common to both 2D and 3D ray Tracing.
 */
#include "common_helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mtwister.h"
#include <mex.h>

/*
 * Linearise [row][column] coordinates in an array of coordinates
 * such as a 3-column n-row matrix of vertices of a Surface3D
 *
 * INPUTS:
 *  row - row indez
 *  col - column index
 *  ind - pointer for the output index to be added to
 */
void lin(int row, int col, int * const ind) {
    *ind = 3*row + col;
}

/*
 * Prints a 2 or 3 element 1D array to the terminal.
 *
 * INPUTS:
 *  vect - a double pointer to either a 2 or 3 element vector
 *  dim  - int, the dimension of the vector, 2 or 3
 */
void print1D_double(double const * const vect, int dim) {
    mexPrintf("[");
    for (int i = 0; i < dim - 1; i++) {
        mexPrintf("%f,", vect[i]);
    }
    mexPrintf("%f]\n", vect[dim - 1]);
}

/*
 * Prints a 2 or 3 element 1D array to the terminal.
 *
 * INPUTS:
 *  vect - an int pointer to either a 2 or 3 element vector
 *  dim  - int, the dimension of the vector, 2 or 3
 */
void print1D_int(int const * const vect, int dim) {
    mexPrintf("[");
    for (int i = 0; i < dim - 1; i++) {
        mexPrintf("%i,", vect[i]);
    }
    mexPrintf("%i]\n", vect[dim - 1]);
}

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double const matrix[3][3]) {
    int i;
    mexPrintf("{\n");
    for (i = 0; i < 3; i++) {
        mexPrintf("%f, ", matrix[i][0]);
        mexPrintf("%f, ", matrix[i][1]);
        mexPrintf("%f\n", matrix[i][2]);
    }
    mexPrintf("}\n");
}

/* 
 * Generates two gaussian random numbers using the box-muller transform.
 */
void gaussian_random(double mu, double sigma, double Z[2], MTRand * const myrng) {
    double U1, U2;
    
    genRand(myrng, &U1);
    genRand(myrng, &U2);
    
    Z[0] = sqrt(-2*log(U1))*cos(2*M_PI*U2);
    Z[1] = sqrt(-2*log(U1))*sin(2*M_PI*U2);
    
    Z[0] = Z[0]*sigma + mu;
    Z[1] = Z[1]*sigma + mu;
    
    return;
}

/* 
 * Samples from the top tail of a Gaussian distribution. Samples from the
 * distribution and then checks to see if the value is below the cutoff.
 *
 * INPUTS:
 *  mu     - mean of the Gaussian function
 *  sigma  - standard deviation of the Gaussian distribution
 *  cutoff - the cutoff that the random number must be larger than
 *  MTRand - random number generator object
 *  rand1  - pointer to where to store the result
 */
void gaussian_random_tail(double mu, double sigma, double cutoff, MTRand * const myrng,
		double * const rand1) {
    double Z[2] = {0, 0};
    int cnt = 0;

    do {
        if (!(cnt % 2)) {
            gaussian_random(mu, sigma, Z, myrng);
            *rand1 = Z[0];
        } else {
            *rand1 = Z[1];
        }
        cnt++;
    } while (*rand1 < cutoff);
}

/*
 * Create a random int in the desired range, 0 to max-1
 */
void gen_random_int(int max, MTRand * const myrand, int * const randint) {
    double uniform_rand;
    
    genRand(myrand, &uniform_rand);
    uniform_rand = uniform_rand*(double)max;
    *randint = (int)floor(uniform_rand);
}
