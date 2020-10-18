/*
 * Copyright (c) 2018-20, Sam Lambrick, Dan Serment.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * A main MEX function for performing the SHeM Simulation.
 *
 * The calling syntax is:
 *  [counted, killed, numScattersRay]  = tracingMultiGenMex(V, F, N, C, sphere, ...
 *  	plate, mat_names, mat_functions, mat_params, max_scatter, n_rays, ...
 *      source_model, source_parameters);
 * 
 * INPUTS:
 *  V - Vertices of the sample
 *  F - Faces of the sample
 *  N - Normals of the sample
 *  C - compositions (materials) of the sample
 *  sphere - matlab struct array of the parameters for an analytic sphere
 *  plate  - matlab struct array of the parameters for the detectors/pinhole plate
 *  mat_names - names of the materials used
 *  mat_functions - names of the functions used
 *  mat_params - array of parameters for the scattering functions used
 *  max_scatter - maximum allowed sample scattering events
 *  n_rays - number of rays to simulate
 *  source_model - string, the source model to use to generate the rays
 *  source_parameter - array of parameters for the source model
 * 
 * OUTPUTS:
 *  counted - number of detected rays into each detector
 *  killed  - number of rays that had to be stopped
 *  numScattesRay - number of scattering events each detected ray underwent
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
    int const NINPUTS = 13;
    int const NOUTPUTS = 3;

    /* Declare the input variables */
    int ntriag_sample;     /* number of sample triangles */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    char **C;              /* sample material keys, length M */
    Material *M;           /* materials of the sample */
    int n_rays;            /* number of rays */
    int maxScatters;       /* Maximum number of scattering events per ray */
    
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
    SourceParam source;
    int i;

    /* Indexing the surfaces, -1 refers to no surface */
    int sample_index = 0, plate_index = 1, sphere_index = 2;
    
    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /**************************************************************************/

    // TODO: improve the input checking
    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMultiGenMex:nrhs",
        		"%d inputs required for tracingMultiGenMex.", NINPUTS);
    }
    if (nlhs != NOUTPUTS) {
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMultiGenMex:nrhs",
        		"%d outputs required for tracingMultiGenMex.", NOUTPUTS);
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
    
    // TODO: pass through source as a struct?
    get_source(prhs[12], (int)mxGetScalar(prhs[11]), &source);

    /**************************************************************************/
        
    // Seed the random number generator with the current time
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    // Set up the MTwister random number generator
    seedRand(t, &myrng);

    // Put the sample and pinhole plate surface into structs
    // TODO: can we make a sample struct that can be passed from Matlab to C?
    set_up_surface(V, N, F, C, M, num_materials, ntriag_sample, nvert, sample_index, &sample);

    /**************************************************************************/
    
    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[0] = mxCreateNumericMatrix(1, plate.n_detect, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, plate.n_detect*maxScatters, mxINT32_CLASS,
                                    mxREAL);

    /* Pointers to the output matrices so we may change them*/
    cntr_detected = (int32_t*)mxGetData(plhs[0]);
    numScattersRay = (int32_t*)mxGetData(plhs[2]);

    /**************************************************************************/

    // FOR DEBUG

    //print_surface(&sample);


    // END DEBUG

    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    // TODO: move loop into experiments.c file
    for (i = 0; i < n_rays; i++) {
        int detected;

        // TODO: create a source struct
        create_ray(&the_ray, &source, &myrng);

        trace_ray_simple_multi(&the_ray, &killed, cntr_detected, maxScatters,
        		 &sample, &plate, &sphere, &detector, &myrng, &detected);
        
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

    return;
}


