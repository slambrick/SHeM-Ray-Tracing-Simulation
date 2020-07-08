/*
 * Copyright (c) 2018-19, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 *
 * The main MEX function for performing the SHeM Simulation.
 *
 * The calling syntax is:
 *
 *  [] = tracingMex();
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
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    char **C;              /* sample triangle diffuse level, length M */
    Material *M;           /* sample scattering parameters */
    int nrays;             /* number of rays */
    int nvert;             /* number of vertices in the sample */
    int ntriag_sample;     /* number of sample triangles */
    int maxScatters;       /* Maximum number of scattering events per ray */
    double *VS;            /* pinhole plate triangle vertices */
    double *FS;            /* pinhole plate triangle faces */
    double *NS;            /* pinhole plate triangle normals */
    double *CS;            /* pinhole plate triangle diffuse level*/
    double *PS;            /* pinhole plate scattering parameters */
    int ntriag_plate;      /* number of pinhole plate triangles */
    double *backWall;
    int make_sphere;       /* Should the analytic sphere be added to the model */
    double *sphere_c;
    double sphere_r;       /* Radius of the analytic sphere if it to be made */
    double sphere_diffuse; /* The scattering off of the analytic sphere */
    double sphere_parameters; /* Scattering distribution parameters */
    int source_model;
    double *source_parameters;

    /* Declare the output variables */
    int cntr_detected;       /* The number of detected rays */
    int killed;              /* The number of killed rays */
    int *numScattersRay;     /* The number of sample scatters that each
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

    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /*******************************************************************************/

    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Twenty inputs required for tracingMex.");
    }
    if (nlhs != NOUTPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Three outpus required for tracingMex.");
    }

    /**************************************************************************/

    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    V = mxGetPr(prhs[0]);
    F = mxGetInt32s(prhs[1]);
    N = mxGetPr(prhs[2]);
    maxScatters = (int)mxGetScalar(prhs[5]); /* mxGetScalar gives a double */
    VS = mxGetPr(prhs[6]);
    FS = mxGetPr(prhs[7]);
    NS = mxGetPr(prhs[8]);
    CS = mxGetPr(prhs[9]);
    PS = mxGetPr(prhs[10]);
    nrays = (int)mxGetScalar(prhs[13]);
    source_model = (int)mxGetScalar(prhs[14]);
    nvert = mxGetN(prhs[0]);
    ntriag_sample = mxGetN(prhs[1]);
    ntriag_plate = mxGetN(prhs[7]);

    // read in the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[3], C);
    
    // get the sphere from struct
    sphere = get_sphere(prhs[11], sphere_index);
    
    // materials
    num_materials = mxGetN(prhs[6]);
    M = mxCalloc(num_materials, sizeof(Material));
    get_materials_array(prhs[6], prhs[7], prhs[8], M);
    
    // extract plate properties from theplate cell array containing plate options
    plate = get_plate(prhs[12], plate_index);
    
    // TODO: pass through source as a struct?
    get_source(prhs[15], &pinhole_r, pinhole_c, &src_theta_max, &src_init_angle,
               &src_sigma);
    
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
    srand(t);
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);

    /* Put the sample and pinhole plate surface into structs */
    // TODO: can we make a sample struct that can be passed from Matlab to C?
    sample = set_up_surface(V, N, F, C, M, num_materials, ntriag_sample, nvert,
                            sample_index);
    sample = set_up_surface(VS, NS, FS, CS, MS, num_materials, ntriag_sample, nvert,
                            sample_index);
    
    /* Put information on the analytic sphere into a struct */
    the_sphere = set_up_sphere(make_sphere, sphere_c,
        sphere_r, sphere_diffuse, sphere_parameters, sphere_index);

    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[2] = mxCreateNumericMatrix(1, maxScatters, mxINT32_CLASS, mxREAL);

    /* Pointers to the output matrices so we may change them*/
    numScattersRay  = (int*)mxGetData(plhs[2]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int detected;

        the_ray = create_ray(source_parameters[0], &source_parameters[1],
            source_parameters[4], source_parameters[5], source_model, 
            source_parameters[6], &myrng);

        detected = trace_ray_triag_plate(&the_ray, &killed, &cntr_detected,
            maxScatters, sample, plate, the_sphere, backWall, &myrng);

        /*
         * Add the number of scattering events the ray has undergon to the
         * histogram. But only if it is detected.
         */
        if (detected)
            numScattersRay[the_ray.nScatters - 1]++;
    }
    /**************************************************************************/

    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(cntr_detected);
    plhs[1] = mxCreateDoubleScalar(killed);

    /* Free the space used by the random number generator */
    //gsl_rng_free(my_rng);

    return;
}
