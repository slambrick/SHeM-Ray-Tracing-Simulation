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
 *  [cntr, killed, final_pos, final_dir, numScattersRay] = ...
 *      tracingMex(init_pos, init_dir, V, F, N, C, maxScatter, VS, FS, ...
 *                 NS, CS, backWall, scan_pos_x, scan_pos_z, make_sphere, ... 
 *                 dist_to_sample, sphere_r, diffuse, first_plate);
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
 *   VS          - lists vertices that make up the pinhole plate
 *   FS          - lists faces that make up the pinhole plate
 *   NS          - lists normals thaat make up the pinhole plate
 *   CS          - indices of the scattering off of each triangle in the 
 *                 pinhole plate
 *   backWall    - the size of the pinhole plate, (y,x,z), don't ask why...
 *   scan_pos_x  - the x scan position (only used when there is a sphere)
 *   scan_pos_z  - the z scan position (only used when there is a sphere)
 *   make_sphere - should an analytic sphere be moddelled
 *   dist_to_sample - the distance between the pinhole plate and the flat 
 *                    surface that the spher sits on (for use when there is a
 *                    sphere)
 *   sphere_r    - the radius of the analytic sphere
 *   diffuse     - index of the scattering off of the analytic sphere
 *   first_plate - should the first scattering event consider the pinhole plate
 *
 *  OUTPUTS:
 *   cntr           - the number of rays that have gone into the detector
 *   killed         - the number of rays that were killed because they reached 
 *                    the max number of scattering events allowed
 *   final_pos      - the final positions of the detected rays
 *   final_dir      - the final directions of the detected rays
 *   numScattersRay - the number of sample scattering events that each ray has 
 *                    undergone
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
//#include <gsl/gsl_rng.h>
#include "small_functions3D.h"
#include "common_helpers.h"
#include "trace_ray.h"
#include "ray_tracing_structs3D.h"
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

/* Function for tracing a single ray */
//int trace_ray(Ray3D *the_ray, int *killed, int *cntr_detected, int maxScatters,
//        Surface3D Sample, Surface3D Plate, AnalytSphere the_sphere, 
//        double backWall[]);

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    
    /* Declare the input variables */
    double *ray_pos;       /* inital ray positions 3xN */
    double *ray_dir;       /* inital ray directions 3xN */
    double *V;             /* sample triangle vertices 3xn */
    double *F;             /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    double *C;             /* sample triangle diffuse level, length M */
    double *P;             /* sample scattering parameters */
    int nrays;             /* number of rays */
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
    //gsl_rng *my_rng;
    int sample_index, sphere_index, plate_index;
    
    /* Declare structs */
    Surface3D Sample;
    Surface3D Plate;
    AnalytSphere the_sphere;
    Rays3D all_rays;
    
    /*******************************************************************************/
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 19) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Nineteen inputs required for tracingMex.");
    }
    if (nlhs != 6) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Six outpus required for tracingMex.");
    }
    
    /* Check the type of the inputs */
    /*if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Input positions must be type double.");
    } if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Input directions must be type double.");
    } if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface vertices must be type double.");
    } if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface face indices must be type double.");
    } if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface normals must be type double.");
    } if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface scattering indices must be type double.");
    } if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]) || 
          mxGetNumberOfElements(prhs[6]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          "Maximum number of sample scatters be real scalar.");
    } if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface vertices must be type double.");
    } if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface face indices must be type double.");
    } if (!mxIsDouble(prhs[9]) || mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface normals must be type double.");
    } if (!mxIsDouble(prhs[10]) || mxIsComplex(prhs[10])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Sample surface scattering indices must be type double.");
    } if (!mxIsDouble(prhs[11]) || mxIsComplex(prhs[11])) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notDouble", 
                          "Information of the detector must be type double.");
    } if (!mxIsDouble(prhs[12]) || mxIsComplex(prhs[12]) || 
          mxGetNumberOfElements(prhs[12]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          "x position must be real scalar.");
    } if (!mxIsDouble(prhs[13]) || mxIsComplex(prhs[13]) || 
          mxGetNumberOfElements(prhs[13]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          "z position must be real scalar.");
    } if (!mxIsDouble(prhs[14]) || mxIsComplex(prhs[14]) || 
          mxGetNumberOfElements(prhs[14]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          " must be real scalar.");
    } if (!mxIsDouble(prhs[15]) || mxIsComplex(prhs[15]) || 
          mxGetNumberOfElements(prhs[15]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          "Distance to the sample must be a real scalar.");
    } if (!mxIsDouble(prhs[16]) || mxIsComplex(prhs[16]) || 
          mxGetNumberOfElements(prhs[16]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          "Sphere radius must be real scalar.");
    } if (!mxIsDouble(prhs[17]) || mxIsComplex(prhs[17]) || 
          mxGetNumberOfElements(prhs[17]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar", 
                          "Scattering idex for the sphere must be real scalar.");
    } if (!mxIsDouble(prhs[18]) || mxIsComplex(prhs[18]) ||
          mxGetNumberOfElements(prhs[18]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:notScalar",
                          "First scattering position must be real scalar.");
    }*/
    
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
    VS = mxGetPr(prhs[8]);
    FS = mxGetPr(prhs[9]);
    NS = mxGetPr(prhs[10]);
    CS = mxGetPr(prhs[11]);
    PS = mxGetPr(prhs[12]);
    backWall = mxGetPr(prhs[13]);
    make_sphere = (int)mxGetScalar(prhs[14]); /* mxGetScalar gives a double */
    sphere_c = mxGetPr(prhs[15]);
    sphere_r = mxGetScalar(prhs[16]);
    sphere_diffuse = mxGetScalar(prhs[17]);
    sphere_parameters = mxGetScalar(prhs[18]);
    nrays = mxGetN(prhs[0]);
    ntriag_sample = mxGetN(prhs[3]);
    ntriag_plate = mxGetN(prhs[9]);
    
    /**************************************************************************/
    
    /* 
     * Number of rays that enter the detector and those that go into the detector
     * after a single scatter 
     */
    cntr_detected = 0;
    
    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;
    
    /* Set up the GSL random number generator */
    //my_rng = setupGSL();
    
    /* See the random number generator with the current time */
    struct timeval tv;
    double t;
    gettimeofday(&tv, 0);
    t =  tv.tv_sec + tv.tv_usec;
    srand(t);
    
    /* Indexing the surfaces, -1 referes to no surface */
    sample_index = 0;
    plate_index = 1;
    sphere_index = 2;
    
    /* Put the rays into a struct */

    /* Put the data into the struct */
    all_rays = compose_rays3D(ray_pos, ray_dir, nrays);
    
    /* Put the sample and pinhole plate surface into structs */
    Sample = set_up_surface(V, N, F, C, P, ntriag_sample, sample_index);
    Plate = set_up_surface(VS, NS, FS, CS, PS, ntriag_plate, plate_index);
    
    /* Put information on the analytic sphere into a struct */
    the_sphere = set_up_sphere(make_sphere, sphere_c,
        sphere_r, sphere_diffuse, sphere_parameters, sphere_index);
    
    plhs[5] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    detected = (int*)mxGetData(plhs[5]);
    
    /**************************************************************************/
    
    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (i = 0; i < all_rays.nrays; i++) {
        detected[i] = trace_ray_triagPlate(&all_rays.rays[i], &killed, &cntr_detected,
            maxScatters, Sample, Plate, the_sphere, backWall);
    }
    
    /**************************************************************************/
    
    /* 
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the 
     * difference in indexing between MATLAB and C.
     */
    plhs[2] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[4] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);

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
    
    /* Free the space used by the random number generator */
    //gsl_rng_free(my_rng);
    
    return;
}
