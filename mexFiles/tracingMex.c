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
    const int NOUTPUTS = 5;

    /* Declare the input variables */
    double *ray_pos;       /* inital ray positions 3xN */
    double *ray_dir;       /* inital ray directions 3xN */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    char **C;              /* sample triangle diffuse level, length M */
    Material *M;           /* sample scattering parameters */
    int nrays;             /* number of rays */
    int nvert_sample;
    int nvert_plate;
    int ntriag_sample;     /* number of sample triangles */
    int maxScatters;       /* Maximum number of scattering events per ray */
    double *VS;            /* pinhole plate triangle vertices */
    int32_t *FS;           /* pinhole plate triangle faces */
    double *NS;            /* pinhole plate triangle normals */
    char **CS;             /* pinhole plate triangle diffuse level*/
    Material *MS;          /* pinhole plate scattering parameters */
    int ntriag_plate;      /* number of pinhole plate triangles */
    double *backWall;

    /* Declare the output variables */
    int cntr_detected;       /* The number of detected rays */
    int killed;              /* The number of killed rays */
    double *final_pos;       /* The final positions of the detected rays */
    double *final_dir;       /* The final directions of the detected rays */
    int *numScattersRay; /* The number of sample scatters that each
                              * ray has undergone */
    int *detected;       /* Logical array, detected? */

    /* Declare other variables */
    int i;
    /* Indexing the surfaces, -1 refers to no surface */
    int sample_index = 0, plate_index = 1, sphere_index = 2;
    int num_materials;

    /* Declare structs */
    Surface3D sample;
    Surface3D plate;
    AnalytSphere the_sphere;
    Rays3D all_rays;

    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /**************************************************************************/

    /* Check for the right number of inputs and outputs */
    if (nrhs != 19) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Nineteen inputs required for tracingMex.");
    }
    if (nlhs != 6) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Six outpus required for tracingMex.");
    }

    /**************************************************************************/

    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    nrays = mxGetN(prhs[0]);
    ray_pos = mxGetPr(prhs[0]);
    ray_dir = mxGetPr(prhs[1]);
    nvert_sample = mxGetN(prhs[2]);
    V = mxGetPr(prhs[2]);
    ntriag_sample = mxGetN(prhs[3]);
    F = mxGetInt32s(prhs[3]); /* Reading in as double not int - cast later in code */
    N = mxGetPr(prhs[4]);

    // read in the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[5], C);

    nvert_plate = mxGetN(prhs[6]);
    VS = mxGetPr(prhs[6]);
    ntriag_plate = mxGetN(prhs[6]);
    FS = mxGetInt32s(prhs[7]);
    NS = mxGetPr(prhs[8]);

    // read in the material keys
    CS = mxCalloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[9], CS);

    // get the sphere from matlab struct array
    the_sphere = get_sphere(prhs[10], sphere_index);

    backWall = mxGetPr(prhs[11]);

    // materials
    int num_materials = mxGetN(prhs[12]);
    M = mxCalloc(num_materials, sizeof(Material));
    get_materials_array(prhs[12], prhs[13], prhs[14], M);

    maxScatters = (int)mxGetScalar(prhs[15]); /* mxGetScalar gives a double */

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

    /* Indexing the surfaces, -1 referes to no surface */
    sample_index = 0;
    plate_index = 1;
    sphere_index = 2;

    /* Put the rays into a struct */
    compose_rays3D(ray_pos, ray_dir, nrays, &all_rays);

    /* Put the sample and pinhole plate surface into structs */
    set_up_surface(V, N, F, C, M, ntriag_sample, nvert_sample, sample_index, &sample);
    set_up_surface(VS, NS, FS, CS, MS, ntriag_plate, nvert_plate, plate_index, &plate);

    /* Put information on the analytic sphere into a struct */
    set_up_sphere(make_sphere, sphere_c, sphere_r, sphere_diffuse, sphere_parameters,
    		sphere_index, &the_sphere);

    plhs[2] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[4] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);

    plhs[5] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    detected = (int*)mxGetData(plhs[5]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */
    given_rays_cad_pinhole(&all_rays, &killed, &cntr_detected, sample, plate,
            the_sphere, backWall, maxScatters, &myrng, numScattersRay);

    /**************************************************************************/

    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */


    /* Pointers to the output matrices so we may change them*/
    final_pos       = mxGetPr(plhs[2]);
    final_dir       = mxGetPr(plhs[3]);
    numScattersRay  = (int*)mxGetData(plhs[4]);

    /* Put data into the output matrices */
    get_positions(&all_rays, final_pos);
    get_directions(&all_rays, final_dir);
    get_scatters(&all_rays, numScattersRay);

    /* Free the allocated memory associated with the rays */
    clean_up_rays(all_rays);

    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(cntr_detected);
    plhs[1] = mxCreateDoubleScalar(killed);

    return;
}
