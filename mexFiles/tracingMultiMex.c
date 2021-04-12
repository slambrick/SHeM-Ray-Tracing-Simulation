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
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mtwister.h"
#include "extract_inputs.h"
#include "atom_ray_tracing3D.h"

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    /* Expected number of inputs and outputs */
    const int NINPUTS = 13;
    const int NOUTPUTS = 5;
    
    /* Declare the input variables */
    double * ray_pos;        /* inital ray positions 3xN */
    double * ray_dir;        /* inital ray directions 3xN */
    double * V;              /* sample triangle vertices 3xn */
    int32_t * F;             /* sample triangle faces 3xM */
    double * N;              /* sample triangle normals 3xM */
    double * B;
    char ** C;               /* sample material keys, length M */
    Material * M;            /* materials of the sample */
    int nrays;               /* number of rays */
    int nvert;               /* number of vertices in the sample */
    int ntriag_sample;       /* number of sample triangles */
    int maxScatters;         /* Maximum number of scattering events per ray */
    
    /* Declare the output variables */
    int32_t * cntr_detected; /* The number of detected rays */
    int killed;              /* The number of killed rays */
    int32_t * numScattersRay;/* The number of sample scatters that each
                              * ray has undergone */
    int * detected;          /* Logical array, detected? */
    int * which_detector;    /* Which detector was the ray detected in */

    /* Declare other variables */
    int sample_index = 0, plate_index = 1, sphere_index = 2;

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
    ray_pos = mxGetDoubles(prhs[0]);
    ray_dir = mxGetDoubles(prhs[1]);
    nvert = mxGetN(prhs[2]);
    V = mxGetDoubles(prhs[2]);
    ntriag_sample = mxGetN(prhs[3]);
    F = mxGetInt32s(prhs[3]);
    N = mxGetDoubles(prhs[4]);
    B = mxGetDoubles(prhs[5]);
    
    // read in the material keys
    C = calloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[6], C);

    // get the sphere from struct
    sphere = get_sphere(prhs[7], sphere_index);

    // extract plate properties from theplate cell array containing plate options
    plate = get_plate(prhs[8], plate_index);

    // materials
    int num_materials = mxGetN(prhs[9]);
    M = calloc(num_materials, sizeof(Material));
    get_materials_array(prhs[9], prhs[10], prhs[11], M);
    
    // simulation parameters
    maxScatters = (int)mxGetScalar(prhs[12]); /* mxGetScalar gives a double */
    
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
    set_up_surface(V, N, B, F, C, M, num_materials, ntriag_sample, nvert, sample_index, &sample);

    /* Output matrix for total number of counts */
    plhs[0] = mxCreateNumericMatrix(1, plate.n_detect, mxINT32_CLASS, mxREAL);
    cntr_detected = (int32_t*)mxGetData(plhs[0]);
    
    /* Output matrix for which rays were detected */
    plhs[3] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    detected = (int32_t*)mxGetData(plhs[3]);
    
    /* Output matrix for which detector the rays went into */
    plhs[4] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    which_detector = (int32_t*)mxGetData(plhs[4]);
    
    /**************************************************************************/
    
    /* Main implementation of the ray tracing */
    given_rays_simple_pinhole(&all_rays, &killed, cntr_detected, sample, plate, sphere,
            maxScatters, detected, which_detector, &myrng);
    
    /**************************************************************************/
    
    /* Output the number of rays we forcefully stopped */
    plhs[1] = mxCreateDoubleScalar(killed);
    
    /* Output matrix for the number of scattering events that each ray underwent */
    plhs[2] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    numScattersRay  = (int32_t*)mxGetData(plhs[2]);
    get_scatters(&all_rays, numScattersRay);
    
    /* Free the allocated memory associated with the rays */
    clean_up_rays(all_rays);
    free(C);
    free(M);
    clean_up_surface(&sample);
    
    return;
}
