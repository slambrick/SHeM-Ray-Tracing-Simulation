/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 * 
 * Contains small helper functions common to both 2D and 3D ray Tracing.
 */
#include "mex.h"
#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <sys/time.h>

/* 
 * Sets up the GNU/SL random number generator and returns a pointer to the
 * generator. This must be called before cosineScatter and the generator this 
 * returns should be passed to cosieScatter.
 * 
 * OUTPUTS:
 *  r - a pointer to a GSL random number generator
 */
gsl_rng* setupGSL(void) {
    gsl_rng *r;
    unsigned long int t;
    
    /*t = (unsigned long int)time(NULL);*/
    
    struct timeval tv;
    gettimeofday(&tv, 0);
    t =  tv.tv_sec + tv.tv_usec;
        
    /* We use the standard random number generator algorithm the Mersenne Twister */
    r = gsl_rng_alloc(gsl_rng_mt19937);
    
    gsl_rng_set(r, t);
    
    return r;
}

/* 
 * Prints a 2 or 3 element 1D array to the terminal.
 * 
 * INPUTS:
 *  vect - a double pointer to either a 2 or 3 element vector
 *  dim  - int, the dimension of the vector, 2 or 3
 */
void print1D_double(double *vect, int dim) {
    if (dim == 2)
        mexPrintf("[%f, %f]\n", vect[0], vect[1]);
    else if (dim == 3)
        mexPrintf("[%f, %f, %f]\n", vect[0], vect[1], vect[2]);
    else
        mexPrintf("dim into print1D must be 2 or 3.");
}

/* 
 * Prints a 2 or 3 element 1D array to the terminal.
 * 
 * INPUTS:
 *  vect - an int pointer to either a 2 or 3 element vector
 *  dim  - int, the dimension of the vector, 2 or 3
 */
void print1D_int(double *vect, int dim) {
    if (dim == 2)
        mexPrintf("[%i, %i]\n", vect[0], vect[1]);
    else if (dim == 3)
        mexPrintf("[%i, %i, %i]\n", vect[0], vect[1], vect[2]);
    else
        mexPrintf("dim into print1D must be 2 or 3.");
}

/* Prints out a 3 by 3 double array passed as an argument. */
void print3x3(double matrix[3][3]) {
    int i;
    mexPrintf("{\n");
    for (i = 0; i < 3; i++) {
        mexPrintf("%f, ", matrix[i][0]);
        mexPrintf("%f, ", matrix[i][1]);
        mexPrintf("%f\n", matrix[i][2]);
    }
    mexPrintf("}\n");
}
