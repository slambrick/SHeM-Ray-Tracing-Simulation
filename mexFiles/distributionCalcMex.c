/*
 * Copyright (c) 2019, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * The main MEX function for performing the SHeM Simulation.
 *
 * The calling syntax is:
 *
 * This is a MEX file for MATLAB.
 */

/*
 * tracing_functions.h results also in the including of:
 *  small_functions.h
 *  math.h
 *  stdio.h
 *  stdlib.h
 *  time.h
 */
#include <mex.h>
#include "small_functions3D.h"
#include "common_helpers.h"
#include "trace_ray.h"
#include "ray_tracing_structs3D.h"
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mtwister.h"

/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    /* Declare the input variables */
    double *V;             /* sample triangle vertices 3xn */
    double *F;             /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    double *C;             /* sample triangle diffuse level, length M */
    double *P;             /* sample scattering parameters */
    int nrays;             /* number of rays */
    int ntriag_sample;     /* number of sample triangles */
    int maxScatters;       /* Maximum number of scattering events per ray */
    double *start_pos;
    double *start_dir;

    /* Declare the output variables */
    int killed;              /* The number of killed rays */
    int *numScattersRay; /* The number of sample scatters that each
                              * ray has undergone */
    double *final_pos;       /* The final positions of the rays */
    double *final_dir;       /* The final directions of the rays */

    /* Declare other variables */
    int i;
    //gsl_rng *my_rng;
    int sample_index;

    /* Declare structs */
    Surface3D Sample;
    AnalytSphere the_sphere;

    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /*******************************************************************************/

    /* Check for the right number of inputs and outputs */
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Twenty inputs required for tracingMex.");
    }
    if (nlhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Three outpus required for tracingMex.");
    }

    /**************************************************************************/

    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    V = mxGetPr(prhs[0]);
    F = mxGetPr(prhs[1]); /* Reading in as double not int - cast later in code */
    N = mxGetPr(prhs[2]);
    C = mxGetPr(prhs[3]);
    P = mxGetPr(prhs[4]);
    maxScatters = (int)mxGetScalar(prhs[5]); /* mxGetScalar gives a double */
    nrays = (int)mxGetScalar(prhs[6]);
    start_pos = mxGetPr(prhs[7]);
    start_dir = mxGetPr(prhs[8]);
    ntriag_sample = mxGetN(prhs[1]);

    /**************************************************************************/

    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;

    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    srand(t);
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);

    /* Indexing the surfaces, -1 referes to no surface */
    sample_index = 0;

    /* Put the sample into a struct */
    Sample = set_up_surface(V, N, F, C, P, ntriag_sample, sample_index);

    /* Define sphere not to exist */
    the_sphere = set_up_sphere(0, start_pos, 1, 1, 1, 2);

    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[1] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3, nrays, mxREAL);

    /* Pointers to the output matrices so we may change them*/
    numScattersRay = (int*)mxGetData(plhs[1]);
    final_pos = mxGetPr(plhs[2]);
    final_dir = mxGetPr(plhs[3]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int j;

        /* Manually initiate the ray as disired */
        the_ray.nScatters = 0;
        the_ray.on_element = -1;
        the_ray.on_surface = -1;
        for (j = 0; j < 3; j++) {
            the_ray.position[j] = start_pos[j];
            the_ray.direction[j] = start_dir[j];
        }

        trace_ray_just_sample(&the_ray, &killed, maxScatters, Sample, the_sphere, &myrng);

        /* Update final position an directions of ray */
        for (j = 0; j < 3; j++) {
            int n;
            n = 3*i + j;
            final_pos[n] = the_ray.position[j];
            final_dir[n] = the_ray.direction[j];
        }
        numScattersRay[i] = the_ray.nScatters;
    }
    /**************************************************************************/

    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(killed);

    /* Free the space used by the random number generator */
    //gsl_rng_free(my_rng);

    return;
}
