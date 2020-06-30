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
 *  [cntr, killed, numScattersRay, detected, which_detector]  = ...
 *      tracingMultiMex(init_pos, init_dir, V, F, N, C, P, maxScatter, ...
 *                 make_sphere, sphere_c, sphere_r, sphere_diffuse, sphere_parameters, ...
 *                 plate_represent, n_detector, circle_plate_r, aperture_axes, aperture_c);
 *  
 *  INPUTS:
 *   init_pos    - the initial positions of the rays
 *   init_dir    - the initial directions of the rays
 *   V           - a list of the locations of the vertices in the surface
 *   F           - lists which vertices make up the triangles in the surface
 *   N           - lists the normals to the triangles in the surface
 *   C           - indices of the scattering off of each triangle in the 
 *                 sample surface
 *   maxScatter  - the maximum number of scatters off of the sample that rays
 *                 are allowed to undergo
 *   make_sphere - should an analytic sphere be moddelled
 *   sphere_c    - the centre coordinate of the analytic sphere
 *   sphere_r    - the radius of the analytic sphere
 *   sphere_diffuse
 *   sphere_parameters
 *   plate_represent
 *   n_detector
 *   circle_plate_r
 *   aperture_axes
 *   aperture_c
 *
 *  OUTPUTS:
 *   cntr           - the number of rays that have gone into the detector
 *   killed         - the number of rays that were killed because they reached 
 *                    the max number of scattering events allowed
 *   numScattersRay - the number of sample scattering events that each ray has 
 *                    undergone
 *   detected       - array, 1/0, is each ray detected
 *   which_detector - array, if the ray was detected which detector did it go 
 *                    into
 * 
 * This is a MEX file for MATLAB.
 */

#include "mex.h"
#include "small_functions3D.h"
#include "common_helpers.h"
#include "trace_ray.h"
#include "ray_tracing_structs3D.h"
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include "mtwister.h"

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    
    /* Declare the input variables */
    double *ray_pos;        /* inital ray positions 3xN */
    double *ray_dir;        /* inital ray directions 3xN */
    double *V;              /* sample triangle vertices 3xn */
    double *F;              /* sample triangle faces 3xM */
    double *N;              /* sample triangle normals 3xM */
    double *C;              /* sample triangle diffuse level, length M */
    double *P;              /* sample scattering parameters */
    int nrays;              /* number of rays */
    int ntriag_sample;      /* number of sample triangles */
    int maxScatters;        /* Maximum number of scattering events per ray */
    int make_sphere;        /* Should the analytic sphere be added to the model */
    double *sphere_c;
    double sphere_r;        /* Radius of the analytic sphere if it to be made */
    double sphere_diffuse;  /* The scattering off of the analytic sphere */
    double sphere_parameters; /* Scattering distribution parameters */
    int plate_represent;    /* Should rays scatter off the pinhole plate */
    int n_detector;         /* How many detectors */
    double circle_plate_r;  /* The radius of the pinhole plate */
    double *aperture_axes;  /* Array of the detector aperture axes */
    double *aperture_c;     /* Array of the centres of the detector apertures */
    
    /* Declare the output variables */
    int *cntr_detected;     /* The number of detected rays */
    int killed;             /* The number of killed rays */
    int *numScattersRay;    /* The number of sample scatters that each
                             * ray has undergone */
    int *detected;          /* Logical array, detected? */
    int *which_detector;    /* Which detector was the ray detected in */

    /* Declare other variables */
    int i;
    int sample_index, sphere_index, plate_index;
    int detector;

    /* Declare structs */
    Surface3D Sample;
    NBackWall Plate;
    AnalytSphere the_sphere;
    Rays3D all_rays;
    
    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /*******************************************************************************/
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 18) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Nineteen inputs required for tracingMex.");
    }
    if (nlhs != 5) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Five outpus required for tracingMex.");
    }
    
    /**************************************************************************/
    
    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    ray_pos = mxGetPr(prhs[0]);
    ray_dir = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[2]);
    F = mxGetPr(prhs[3]); /* Reading in as double not int - cast later in code */
    N = mxGetPr(prhs[4]);
    C = mxGetPr(prhs[5]);
    P = mxGetPr(prhs[6]);
    maxScatters = (int)mxGetScalar(prhs[7]); /* mxGetScalar gives a double */
    make_sphere = (int)mxGetScalar(prhs[8]); /* mxGetScalar gives a double */
    sphere_c = mxGetPr(prhs[9]);
    sphere_r = mxGetScalar(prhs[10]);
    sphere_diffuse = mxGetScalar(prhs[11]);
    sphere_parameters = mxGetScalar(prhs[12]);
    plate_represent = (int)mxGetScalar(prhs[13]);
    n_detector = (int)mxGetScalar(prhs[14]);
    circle_plate_r = mxGetScalar(prhs[15]);
    aperture_axes = mxGetPr(prhs[16]);
    aperture_c = mxGetPr(prhs[17]);
    nrays = mxGetN(prhs[0]);
    ntriag_sample = mxGetN(prhs[3]);
    
    /**************************************************************************/

    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    srand(t);
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);
    
    /* Indexing the surfaces, -1 referes to no surface */
    sample_index = 0;
    plate_index = 1;
    sphere_index = 2;
    
    /* Put the rays into a struct */
    all_rays = compose_rays3D(ray_pos, ray_dir, nrays);
    
    /* Put the sample and pinhole plate surface into structs */
    Sample = set_up_surface(V, N, F, C, P, ntriag_sample, sample_index);
    Plate.aperture_c = aperture_c;
    Plate.aperture_axes = aperture_axes;
    Plate.circle_plate_r = circle_plate_r;
    Plate.composition = 1;
    Plate.scattering_parameters = 0;
    Plate.plate_represent = plate_represent;
    Plate.surf_index = plate_index;
    Plate.n_detect = n_detector;
    
    /* Put information on the analytic sphere into a struct */
    the_sphere = set_up_sphere(make_sphere, sphere_c,
        sphere_r, sphere_diffuse, sphere_parameters, sphere_index);
    
    /* Output matrix for total number of counts */
    plhs[0] = mxCreateNumericMatrix(1, n_detector, mxINT32_CLASS, mxREAL);
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
    for (i = 0; i < all_rays.nrays; i++) {
        detected[i] = trace_ray_simpleMulti(&all_rays.rays[i], &killed, cntr_detected,
            maxScatters, Sample, Plate, the_sphere, &detector, &myrng);
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
    
    return;
}
