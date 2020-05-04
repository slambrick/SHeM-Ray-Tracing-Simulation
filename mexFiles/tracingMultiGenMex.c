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

#include "mex.h"
#include "trace_ray.h"
#include "small_functions3D.h"
#include "common_helpers.h"
#include "ray_tracing_structs3D.h"
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    
    /* Declare the input variables */
    double *V;             /* sample triangle vertices 3xn */
    double *F;             /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    double *C;             /* sample triangle diffuse level, length M */
    double *P;             /* sample scattering parameters */
    int nrays;             /* number of rays */
    int ntriag_sample;     /* number of sample triangles */
    int maxScatters;       /* Maximum number of scattering events per ray */
    int make_sphere;       /* Should the analytic sphere be added to the model */
    double *sphere_c;     /* The pinhole-plate sample distance (for use with 
                            * the analytic sphere) */
    double sphere_r;       /* Radius of the analytic sphere if it to be made */
    double sphere_diffuse; /* The scattering off of the analytic sphere */
    double sphere_parameters; /* Scattering distribution parameters */
    int plate_represent;
    int n_detector;
    double circle_plate_r;
    double *aperture_axes;
    double *aperture_c;
    int source_model;
    double *source_parameters;
    
    /* Declare the output variables */
    int *cntr_detected;       /* The number of detected rays */
    int killed;              /* The number of killed rays */
    int *numScattersRay; /* The number of sample scatters that each
                              * ray has undergone */
    
    /* Declare other variables */
    int i;
    int sample_index, sphere_index, plate_index;
    int detector;
    
    /* Declare structs */
    Surface3D Sample;
    NBackWall Plate;
    AnalytSphere the_sphere;
    
    /*******************************************************************************/
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 19) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Nineteen inputs required for tracingMex.");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Six outpus required for tracingMex.");
    }
    
    /**************************************************************************/
    
    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    V = mxGetPr(prhs[0]);
    F = mxGetPr(prhs[1]); /* Reading in as double not int - cast later in code */
    N = mxGetPr(prhs[2]);
    C = mxGetPr(prhs[3]);
    P = mxGetPr(prhs[4]);
    maxScatters = (int)mxGetScalar(prhs[5]); /* mxGetScalar gives a double */
    make_sphere = (int)mxGetScalar(prhs[6]); /* mxGetScalar gives a double */
    sphere_c = mxGetPr(prhs[7]);
    sphere_r = mxGetScalar(prhs[8]);
    sphere_diffuse = mxGetScalar(prhs[9]);
    sphere_parameters = mxGetScalar(prhs[10]);
    plate_represent = (int)mxGetScalar(prhs[11]);
    n_detector = (int)mxGetScalar(prhs[12]);
    circle_plate_r = mxGetScalar(prhs[13]);
    aperture_axes = mxGetPr(prhs[14]);
    aperture_c = mxGetPr(prhs[15]);
    nrays = (int)mxGetScalar(prhs[16]);
    source_model = (int)mxGetScalar(prhs[17]);
    source_parameters = mxGetPr(prhs[18]);
    ntriag_sample = mxGetN(prhs[1]);
    
    
    /**************************************************************************/
    
    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;
    
    /* Set up the GSL random number generator */
    /*my_rng = setupGSL();*/
    
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
        
    /* 
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the 
     * difference in indexing between MATLAB and C.
     */
    plhs[0] = mxCreateNumericMatrix(1, n_detector, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, n_detector*maxScatters, mxINT32_CLASS, mxREAL);
    
    /* Pointers to the output matrices so we may change them*/
    cntr_detected = (int*)mxGetData(plhs[0]);
    numScattersRay = (int*)mxGetData(plhs[2]);
    
    /**************************************************************************/
    
    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int detected;
        
        the_ray = create_ray_source(source_parameters[0], &source_parameters[1], 
            source_parameters[4], source_parameters[5], source_model, 
            source_parameters[6]);
        
        detected = trace_ray_simpleMulti(&the_ray, &killed, cntr_detected,
            maxScatters, Sample, Plate, the_sphere, &detector);
                
        /* 
         * Add the number of scattering events the ray has undergon to the
         * histogram. But only if it is detected.
         */
        if (detected) {
            int ind;
            ind = (detector - 1)*maxScatters + (the_ray.nScatters - 1);
            numScattersRay[ind]++;
        }
    }
    
    /**************************************************************************/
    
    /* Output number of rays went into the detector */
    plhs[1] = mxCreateDoubleScalar(killed);
    
    return;
}


