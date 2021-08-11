/*
 * Copyright (c) 2018-20, Sam Lambrick.
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
    const int NINPUTS = 20;
    const int NOUTPUTS = 6;

    /* Declare the input variables */
    double *ray_pos;       /* inital ray positions 3xN */
    double *ray_dir;       /* inital ray directions 3xN */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    double *B;
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
    double *BS;
    char **CS;             /* pinhole plate triangle diffuse level*/
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

    /* Indexing the surfaces, -1 refers to no surface */
    int sample_index = 0, plate_index = 1, circle_index = 2, sphere_index = 3;

    /* Declare structs */
    Surface3D sample;
    Surface3D plate;
    AnalytSphere * spheres;
    Circle the_circle;
    Rays3D all_rays;
    Sample overall_sample;
    int n_sphere;

    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /**************************************************************************/

    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "Sizteen inputs required for tracingMex.");
    }
    if (nlhs != NOUTPUTS) {
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
    B = mxGetPr(prhs[5]);

    // read in the material keys
    C = calloc(ntriag_sample, sizeof(char*));
    get_string_cell_arr(prhs[6], C);

    nvert_plate = mxGetN(prhs[7]);
    VS = mxGetPr(prhs[7]);
    ntriag_plate = mxGetN(prhs[8]);
    FS = mxGetInt32s(prhs[8]);
    NS = mxGetPr(prhs[9]);
    BS = mxGetPr(prhs[10]);

    // read in the material keys
    CS = calloc(ntriag_plate, sizeof(char*));
    get_string_cell_arr(prhs[11], CS);

    // get the sphere from matlab struct array
    // TODO: actual make this work
    n_sphere = (int)mxGetScalar(prhs[12]);
    spheres = (AnalytSphere *)malloc(n_sphere*sizeof(AnalytSphere));
    get_spheres(n_sphere, prhs[13], sphere_index, spheres);
    //the_sphere = get_sphere(prhs[12], sphere_index);
    the_circle = get_circle(prhs[14], circle_index);

    backWall = mxGetPr(prhs[15]);

    // materials
    int num_materials = mxGetN(prhs[16]);
    M = calloc(num_materials, sizeof(Material));
    get_materials_array(prhs[16], prhs[17], prhs[18], M);

    maxScatters = (int)mxGetScalar(prhs[19]); /* mxGetScalar gives a double */

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

    /* Put the rays into a struct */
    compose_rays3D(ray_pos, ray_dir, nrays, &all_rays);

    /* Put the sample and pinhole plate surface into structs */
    set_up_surface(V, N, B, F, C, M, num_materials, ntriag_sample, nvert_sample, sample_index, &sample);
    set_up_surface(VS, NS, BS, FS, CS, M, num_materials, ntriag_plate, nvert_plate, plate_index, &plate);

    /* Put all the sample structs together in one struct */
    overall_sample.the_sphere = spheres;
    overall_sample.the_circle = &the_circle;
    overall_sample.triag_sample = &sample;
    overall_sample.n_sphere = n_sphere;
    
    plhs[2] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    final_pos = (double *)mxGetData(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    final_dir = (double *)mxGetPr(plhs[3]);
    plhs[4] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    numScattersRay  = (int32_t *)mxGetData(plhs[4]);

    plhs[5] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    detected = (int32_t *)mxGetData(plhs[5]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */
    given_rays_cad_pinhole(&all_rays, &killed, &cntr_detected, overall_sample, plate,
            backWall, maxScatters, detected, &myrng);

    /**************************************************************************/

    /* Put data into the output matrices */
    get_positions(&all_rays, final_pos);
    get_directions(&all_rays, final_dir);
    get_scatters(&all_rays, numScattersRay);

    /* Free space */
    free(C);
    free(CS);
    free(M);
    clean_up_surface(&sample);
    clean_up_surface(&plate);
    clean_up_rays(all_rays);

    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(cntr_detected);
    plhs[1] = mxCreateDoubleScalar(killed);

    return;
}
