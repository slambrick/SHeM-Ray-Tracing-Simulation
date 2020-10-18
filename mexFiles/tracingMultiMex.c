/*
 * Copyright (c) 2020, Sam Lambrick.
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

#include <mex.h>
#include <matrix.h>
#include <stdint-gcc.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mtwister.h"
#include "trace_ray.h"
#include "extract_inputs.h"
#include "common_helpers.h"
#include "ray_tracing_structs3D.h"

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    /* Expected number of inputs and outputs */
    const int NINPUTS = 12;
    const int NOUTPUTS = 5;
    
    /* Declare the input variables */
    double *ray_pos;        /* inital ray positions 3xN */
    double *ray_dir;        /* inital ray directions 3xN */
    double *V;              /* sample triangle vertices 3xn */
    int32_t *F;              /* sample triangle faces 3xM */
    double *N;              /* sample triangle normals 3xM */
    char **C;              /* sample triangle diffuse level, length M */
    Material *M;              /* sample scattering parameters */
    int nrays;              /* number of rays */
    int nvert;              /* number of vertices in the sample */
    int ntriag_sample;      /* number of sample triangles */
    int maxScatters;        /* Maximum number of scattering events per ray */
    
    /* Declare the output variables */
    int *cntr_detected;     /* The number of detected rays */
    int killed;             /* The number of killed rays */
    int *numScattersRay;    /* The number of sample scatters that each
                             * ray has undergone */
    int *detected;          /* Logical array, detected? */
    int *which_detector;    /* Which detector was the ray detected in */

    /* Declare other variables */
    int i;
    int sample_index = 0, plate_index = 1, sphere_index = 2;
    int detector;

    /* Declare structs */
    Surface3D sample;
    NBackWall plate;
    AnalytSphere sphere;
    Rays3D all_rays;
    
    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /*******************************************************************************/
    
    // TODO: improve the input checking
    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMultiMex:nrhs",
        		"%d inputs required for tracingMultiMex.", NINPUTS);
    }
    if (nlhs != NOUTPUTS) {
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMultiMex:nrhs",
        		"%d outputs required for tracingMultiMex.", NOUTPUTS);
    }
    
    /**************************************************************************/
    
    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    nrays = mxGetN(prhs[0]);
    ray_pos = mxGetPr(prhs[0]);
    ray_dir = mxGetPr(prhs[1]);
    nvert = mxGetN(prhs[2]);
    V = mxGetPr(prhs[2]);
    ntriag_sample = mxGetN(prhs[3]);
    F = mxGetInt32s(prhs[3]);
    N = mxGetPr(prhs[4]);
    
    // read in the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[5], C);

    // get the sphere from struct
    sphere = get_sphere(prhs[6], sphere_index);

    // extract plate properties from theplate cell array containing plate options
    plate = get_plate(prhs[7], plate_index);

    // materials
    int num_materials = mxGetN(prhs[8]);
    M = mxCalloc(num_materials, sizeof(Material));
    get_materials_array(prhs[8], prhs[9], prhs[10], M);
    
    // simulation parameters
    maxScatters = (int)mxGetScalar(prhs[11]); /* mxGetScalar gives a double */
    nrays = (int)mxGetScalar(prhs[12]);
    
    /**************************************************************************/

    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    /* Set up the MTwister random number generator */
    seedRand(t, &myrng);

    /* Put the rays into a struct */
    compose_rays3D(ray_pos, ray_dir, nrays, &all_rays);
    
    /* Put the sample and pinhole plate surface into structs */
    set_up_surface(V, N, F, C, M, num_materials, ntriag_sample, nvert, sample_index, &sample);
    
    /* Output matrix for total number of counts */
    plhs[0] = mxCreateNumericMatrix(1, plate.n_detect, mxINT32_CLASS, mxREAL);
    cntr_detected = (int*)mxGetData(plhs[0]);
    
    /* Output matrix for which rays were detected */
    plhs[3] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    detected = (int*)mxGetData(plhs[3]);
    
    /* Output matrix for which detector the rays went into */
    plhs[4] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    which_detector = (int*)mxGetData(plhs[4]);
    
    /**************************************************************************/
    
    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    // TODO: move loop into experiments.c file
    for (i = 0; i < all_rays.nrays; i++) {
        trace_ray_simple_multi(&all_rays.rays[i], &killed, cntr_detected, maxScatters,
        		&sample, &plate, &sphere, &detector, &myrng, &detected[i]);
        which_detector[i] = detector;
    }
    
    /**************************************************************************/
    
    /* Output the number of rays we forcefully stopped */
    plhs[1] = mxCreateDoubleScalar(killed);
    
    /* Output matrix for the number of scattering events that each ray underwemt */
    plhs[2] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    numScattersRay  = (int*)mxGetData(plhs[2]);
    get_scatters(&all_rays, numScattersRay);
    
    /* Free the allocated memory associated with the rays */
    clean_up_rays(all_rays);
    mxFree(C);
    mxFree(M);
    clean_up_surface(&sample);
    
    return;
}
