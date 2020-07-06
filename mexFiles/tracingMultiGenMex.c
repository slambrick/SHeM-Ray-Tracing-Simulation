/*
 * Copyright (c) 2018-20, Sam Lambrick, Dan Serment.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * The main MEX function for performing the SHeM Simulation.
 *
 * The calling syntax is:
 *
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
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Expected number of inputs and outputs */
    const int NINPUTS = 13;
    const int NOUTPUTS = 3;

    /* Declare the input variables */
    int ntriag_sample;     /* number of sample triangles */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    char **C;              /* sample material keys, length M */
    Material *M;           /* materials of the sample */
    int n_rays;            /* number of rays */
    int maxScatters;       /* Maximum number of scattering events per ray */
    int source_model;
    double pinhole_r;
    double src_theta_max;
    double src_init_angle;
    double src_sigma;
    double pinhole_c[3];
    
    /* Declare the output variables */
    int32_t * cntr_detected;       /* The number of detected rays */
    int killed = 0;                /* The number of killed rays */
    int32_t * numScattersRay;      /* The number of sample scatters that each
                                    * ray has undergone */

    /* Declare other variables */
    int detector = 0;
    int nvert;
    Surface3D sample;
    NBackWall plate;
    AnalytSphere sphere;
    Ray3D the_ray;

    /* Indexing the surfaces, -1 refers to no surface */
    int sample_index = 0, plate_index = 1, sphere_index = 2;
    
    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /**************************************************************************/

    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "%d inputs required for tracingMultiGenMex.", NINPUTS);
    }
    if (nlhs != NOUTPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "%d outpus required for tracingMultiGenMex.", NOUTPUTS);
    }

    /**************************************************************************/
    
    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    nvert = mxGetN(prhs[0]);
    V = mxGetDoubles(prhs[0]);

    ntriag_sample = mxGetN(prhs[1]);
    F = mxGetInt32s(prhs[1]);
    N = mxGetDoubles(prhs[2]);

    // read in the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[3], C);

    // get the sphere from struct
    sphere = get_sphere(prhs[4], sphere_index);

    // extract plate properties from thePlate cell array containing plate options
    plate = get_plate(prhs[5], plate_index);

    // materials
    int num_materials = mxGetN(prhs[6]);
    M = mxCalloc(num_materials, sizeof(Material));
    get_materials_array(prhs[6], prhs[7], prhs[8], M);
    
    // simulation parameters
    maxScatters = (int)mxGetScalar(prhs[9]);
    n_rays = (int)mxGetScalar(prhs[10]);
    source_model = (int)mxGetScalar(prhs[11]);
    
    // TODO: pass through source as a struct?
    get_source(prhs[12], &pinhole_r, pinhole_c, &src_theta_max, &src_init_angle, &src_sigma);

    /**************************************************************************/
        
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);

    /* Put the sample and pinhole plate surface into structs */
    // TODO: can we make a sample struct that can be passed from Matlab to C?
    sample = set_up_surface(V, N, F, C, M, num_materials, ntriag_sample, nvert, sample_index);

    /**************************************************************************/
    
    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[0] = mxCreateNumericMatrix(1, plate.n_detect, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, plate.n_detect*maxScatters, mxINT32_CLASS, mxREAL);

    /* Pointers to the output matrices so we may change them*/
    cntr_detected = (int*)mxGetData(plhs[0]);
    numScattersRay = (int*)mxGetData(plhs[2]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (int i = 0; i < n_rays; i++) {
        int detected;

        // TODO: create a source struct
        create_ray(&the_ray, pinhole_r, pinhole_c, src_theta_max,
            src_init_angle, source_model, src_sigma, &myrng);

        detected = trace_ray_simple_multi(&the_ray, &killed, cntr_detected,
            maxScatters, sample, plate, sphere, &detector, &myrng);
        
        /*
         * Add the number of scattering events the ray has undergone to the
         * histogram. But only if it is detected.
         */
        if (detected) {
            int ind = (detector - 1)*maxScatters + (the_ray.nScatters - 1);
            numScattersRay[ind]++;
        }
    }

    /**************************************************************************/

    /* Output number of rays went into the detector */
    plhs[1] = mxCreateDoubleScalar(killed);

    /* Free space */
    mxFree(C);
    mxFree(M);
    clean_up_surface(&sample);
}


