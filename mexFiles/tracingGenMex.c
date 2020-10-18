/*
 * Copyright (c) 2018-20, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * A main MEX function for performing the SHeM Simulation.
 *
 * The calling syntax is:
 *
 * [cntr, killed, numScattersRay]  = ...
 *        tracingGenMex(VT, FT, NT, CT, VTS, FTS, NTS, CTS, s, backWall, mat_names, ...
 *                mat_functions, mat_params, max_scatter, beam.n, source_model, source_parameters);
 *
 *  INPUTS:
 *   -
 *
 *  OUTPUTS:
 *   -
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
    const int NINPUTS = 16;
    const int NOUTPUTS = 3;

    /* Declare the input variables */
    double *V;              /* sample triangle vertices 3xn */
    int32_t *F;             /* sample triangle faces 3xM */
    double *N;              /* sample triangle normals 3xM */
    char **C;               /* sample triangle diffuse level, length M */
    Material *M;            /* sample scattering parameters */
    int nrays;              /* number of rays */
    int nvert_sample;       /* number of vertices in the sample */
    int nvert_plate;        /* number of vertices in the plate */
    int ntriag_sample;      /* number of sample triangles */
    int maxScatters;        /* Maximum number of scattering events per ray */
    double *VS;             /* pinhole plate triangle vertices */
    double *FS;             /* pinhole plate triangle faces */
    double *NS;             /* pinhole plate triangle normals */
    char **CS;              /* pinhole plate triangle diffuse level*/
    double *PS;             /* pinhole plate scattering parameters */
    int ntriag_plate;       /* number of pinhole plate triangles */
    double *backWall;

    /* Declare the output variables */
    int cntr_detected;        /* The number of detected rays */
    int killed;               /* The number of killed rays */
    int32_t * numScattersRay; /* The number of sample scatters that each
                              * ray has undergone */

    /* Declare other variables */
    int i;
    /* Indexing the surfaces, -1 refers to no surface */
    int sample_index = 0, plate_index = 1, sphere_index = 2;
    int num_materials;

    /* Declare structs */
    Surface3D sample;
    Surface3D plate;
    AnalytSphere the_sphere;
    SourceParam source;

    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /**************************************************************************/

    // TODO: improve the input checking
    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("AtomRayTracing:tracingGenMex:nrhs",
        		"%d inputs required for tracingGenMex.", NINPUTS);
    }
    if (nlhs != NOUTPUTS) {
        mexErrMsgIdAndTxt("AtomRayTracing:tracingGenMex:nrhs",
        		"%d outputs required for tracingGenMex.", NOUTPUTS);
    }

    /**************************************************************************/

    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    nvert_sample = mxGetN(prhs[0]);
    V = mxGetPr(prhs[0]);
    F = mxGetInt32s(prhs[1]);
    ntriag_sample = mxGetN(prhs[1]);
    N = mxGetPr(prhs[2]);

    // read in the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[3], C);

    nvert_plate = mxGetN(prhs[4]);
    VS = mxGetPr(prhs[4]);
    ntriag_plate = mxGetN(prhs[5]);
    FS = mxGetPr(prhs[5]);
    NS = mxGetPr(prhs[6]);

    // read in the material keys
    CS = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[7], CS);

    // get the sphere from matlab struct array
    the_sphere = get_sphere(prhs[8], sphere_index);

    backWall = mxGetPr(prhs[9]);

    // materials
    int num_materials = mxGetN(prhs[10]);
    M = mxCalloc(num_materials, sizeof(Material));
    get_materials_array(prhs[10], prhs[11], prhs[12], M);

    maxScatters = (int)mxGetScalar(prhs[13]); /* mxGetScalar gives a double */
    nrays = (int)mxGetScalar(prhs[14]);

    // TODO: pass through source as a struct?
    get_source(prhs[16], (int)mxGetScalar(prhs[15]), &source);
    
    /**************************************************************************/

    /*
     * Number of rays that enter the detector and those that go into the detector
     * after a single scatter
     */
    cntr_detected = 0;

    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;

    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;

    /* Set up the MTwister random number generator */
    seedRand(t, &myrng);

    /* Put the sample and pinhole plate surface into structs */
    // TODO: can we make a sample struct that can be passed from Matlab to C?
    set_up_surface(V, N, F, C, M, num_materials, ntriag_sample, nvert_sample, sample_index, &sample);
    set_up_surface(VS, NS, FS, CS, num_materials, ntriag_sample, nvert_plate, plate_index, &plate);

    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[2] = mxCreateNumericMatrix(1, maxScatters, mxINT32_CLASS, mxREAL);

    /* Pointers to the output matrices so we may change them*/
    numScattersRay  = (int32_t*)mxGetData(plhs[2]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    // TODO: move loop into experiments.c file
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int detected;

        create_ray(&the_ray, &source, &myrng);

        trace_ray_triag_plate(&the_ray, &killed, &cntr_detected, maxScatters, &sample,
        		&plate, &the_sphere, backWall, &myrng, &detected);

        if (detected)
            numScattersRay[the_ray.nScatters - 1]++;
    }

    /**************************************************************************/

    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(cntr_detected);
    plhs[1] = mxCreateDoubleScalar(killed);

    /* Free the space used by the random number generator */

    return;
}
