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
#include "mex.h"
#include "trace_ray.h"
#include "extract_inputs.h"
#include "atom_ray_tracing3D.h"
#include "mtwister.h"
#include <stdio.h>
#include <math.h>
#include <stdint-gcc.h>
#include <sys/time.h>

/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    /* Declare the input variables */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    double * B;
    char **C;              /* sample composition */
    Material *M;           /* materials of the sample */
    int nrays;             /* number of rays */
    int nvert;             /* number of sample vertices */
    int ntriag_sample;     /* number of sample triangles */
    int maxScatters;       /* Maximum number of scattering events per ray */
    double *start_pos;
    double *start_dir;
    int n_provided_rays;

    /* Declare the output variables */
    int killed = 0;          /* # rays stopped '.' they scattered too many times */
    int32_t *numScattersRay; /* The number of sample scatters that each
                              * ray has undergone */
    double *final_pos;       /* The final positions of the rays */
    double *final_dir;       /* The final directions of the rays */

    // Declare structs
    Surface3D sample;
    AnalytSphere the_sphere;
    Circle the_circle;
    Sample overall_sample;

    // surface indexing: -1 is no surface, 0 is the sample, etc
    int sample_index = 0;   
    int sphere_index = 1;
    Rays3D all_rays;
    int gen_rays;
    
    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /*******************************************************************************/

    /* Check for the right number of inputs and outputs */
    if (nrhs != 11) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "11 inputs required for distributionCalcMex.");
    }
    if (nlhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "4 outpus required for distributionCalcMex.");
    }

    /**************************************************************************/

    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    V = mxGetPr(prhs[0]);
    F = mxGetInt32s(prhs[1]);
    N = mxGetPr(prhs[2]);
    B = mxGetDoubles(prhs[3]);
    maxScatters = (int)mxGetScalar(prhs[8]); /* mxGetScalar gives a double */
    nrays = (int)mxGetScalar(prhs[9]);
    start_pos = mxGetPr(prhs[10]);
    start_dir = mxGetPr(prhs[11]);
    nvert = mxGetN(prhs[0]);
    ntriag_sample = mxGetN(prhs[1]);
    
    // Get the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[4], C);
    
    // Get the materials
    int num_materials = mxGetN(prhs[4]);
    M = mxCalloc(num_materials, sizeof(Material));
    get_materials_array(prhs[5], prhs[6], prhs[7], M);
    
    /**************************************************************************/
    
    /* Depending on the size of the starting vectors we may or may not generate
     * ray positions ourselves. */
    n_provided_rays = mxGetN(prhs[9]);
    gen_rays = n_provided_rays == 1;
    if (nrays != n_provided_rays) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Provided ray positions is neither 1 or the specified number of rays, for distributionCalcMex.");
    }
    if (!gen_rays) {
        compose_rays3D(start_pos, start_dir, n_provided_rays, &all_rays);
    }
        
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    /* Set up the MTwister random number generator */
    seedRand(t, &myrng);
    
    /* Put the sample into a struct */
    set_up_surface(V, N, B, F, C, M, num_materials, ntriag_sample, nvert, sample_index, &sample);

    /* Define sphere not to exist */
    // TODO: sphere to be passed in as a struct
    set_up_sphere(0, start_pos, 1, M[0], sphere_index, &the_sphere);

    /* Put all the sample structs together in one struct */
    overall_sample.the_sphere = &the_sphere;
    overall_sample.the_circle = &the_circle;
    overall_sample.triag_sample = &sample;
    
    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[1] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3, nrays, mxREAL);

    /* Pointers to the output matrices so we may change them*/
    numScattersRay = (int32_t*)mxGetData(plhs[1]);
    final_pos = mxGetPr(plhs[2]);
    final_dir = mxGetPr(plhs[3]);

    /**************************************************************************/

    /* Loop through all the rays, tracing each one */
    int i, j;
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        
        if (gen_rays) {
            /* Manually initiate the ray as disired */
            the_ray.nScatters = 0;
            the_ray.on_element = -1;
            the_ray.on_surface = -1;
            for (j = 0; j < 3; j++) {
                the_ray.position[j] = start_pos[j];
                the_ray.direction[j] = start_dir[j];
            }
            
            trace_ray_just_sample(&the_ray, &killed, maxScatters, overall_sample, &myrng);
            /* Update final position an directions of ray */
            for (j = 0; j < 3; j++) {
                int n;
                n = 3*i + j;
                final_pos[n] = the_ray.position[j];
                final_dir[n] = the_ray.direction[j];
            }
            numScattersRay[i] = the_ray.nScatters;
        } else {
            trace_ray_just_sample(&all_rays.rays[i], &killed, maxScatters, overall_sample,
                              &myrng);
        }
    }
    
    /**************************************************************************/

    if (!gen_rays) {
        get_scatters(&all_rays, numScattersRay);
        get_positions(&all_rays, final_pos);
        get_directions(&all_rays, final_dir);
        clean_up_rays(all_rays);
    }
    
    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(killed);

    /* Free space */
    mxFree(C);
    mxFree(M);
    clean_up_surface(&sample);
}
