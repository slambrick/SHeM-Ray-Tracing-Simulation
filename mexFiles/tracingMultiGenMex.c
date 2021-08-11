/*
 * Copyright (c) 2018-21, Sam Lambrick, Dan Serment.
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
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mtwister.h"
#include "atom_ray_tracing3D.h"
#include "extract_inputs.h"


/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Expected number of inputs and outputs */
    int const NINPUTS = 16;
    int const NOUTPUTS = 3;

    /* Declare the input variables */
    int ntriag_sample;     /* number of sample triangles */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    double *B;             /* surface reciprocal lattice vectors */
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
    int nvert;
    Surface3D sample;
    NBackWall plate;
    AnalytSphere * spheres;
    Circle the_circle;
    SourceParam source;
    Sample overall_sample;
    int n_sphere;

    /* Indexing the surfaces, -1 refers to no surface */
    /* TODO: make this work a bit better */
    int sample_index = 0, plate_index = 1, circle_index = 2, sphere_index = 3;
    
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
        mexErrMsgIdAndTxt("AtomRayTracing:tracingMultiGenMex:nlhs",
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
    B = mxGetDoubles(prhs[3]);

    // read in the material keys
    C = calloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[4], C);

    // get the sphere from struct
    // TODO: actual make this work
    n_sphere = (int)mxGetScalar(prhs[5]);
    spheres = (AnalytSphere *)malloc(n_sphere*sizeof(AnalytSphere));
    get_spheres(n_sphere, prhs[6], sphere_index, spheres);

    //sphere = get_sphere(prhs[5], sphere_index);
    the_circle = get_circle(prhs[7], circle_index);

    // extract plate properties from thePlate cell array containing plate options
    plate = get_plate(prhs[8], plate_index);
    
    // materials
    int num_materials = mxGetN(prhs[9]);
    M = calloc(num_materials, sizeof(Material));
    get_materials_array(prhs[9], prhs[10], prhs[11], M);
    
    // simulation parameters
    maxScatters = (int)mxGetScalar(prhs[12]);
    n_rays = (int)mxGetScalar(prhs[13]);
    
    // TODO: pass through source as a struct?
    get_source(prhs[15], (int)mxGetScalar(prhs[14]), &source);

    /**************************************************************************/
        
    // Seed the random number generator with the current time
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    // Set up the MTwister random number generator
    seedRand(t, &myrng);

    // Put the sample and pinhole plate surface into structs
    // TODO: can we make a sample struct that can be passed from Matlab to C?
    set_up_surface(V, N, B, F, C, M, num_materials, ntriag_sample, nvert, sample_index, &sample);

    /**************************************************************************/
    
    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[0] = mxCreateNumericMatrix(1, plate.n_detect, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, plate.n_detect*maxScatters, mxINT32_CLASS,
                                    mxREAL);

    /* Put all the sample structs together in one struct */
    overall_sample.the_sphere = spheres;
    overall_sample.the_circle = &the_circle;
    overall_sample.triag_sample = &sample;
    overall_sample.n_sphere = n_sphere;
    
    /* Pointers to the output matrices so we may change them*/
    cntr_detected = (int32_t*)mxGetData(plhs[0]);
    numScattersRay = (int32_t*)mxGetData(plhs[2]);

    /**************************************************************************/

    //print_spheres(spheres, n_sphere);
    /* Main implementation of the ray tracing */
    generating_rays_simple_pinhole(source, n_rays, &killed, cntr_detected,
            maxScatters, overall_sample, plate, &myrng, numScattersRay);

    /**************************************************************************/

    plhs[1] = mxCreateDoubleScalar(killed);

    /* Free space */
    free(C);
    free(M);
    free(spheres);
    clean_up_surface(&sample);

    return;
}


