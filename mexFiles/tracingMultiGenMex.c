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
#include <gsl/gsl_rng.h>
#include "trace_ray.h"
#include "small_functions3D.h"
#include "common_helpers.h"
#include "ray_tracing_structs3D.h"
#include <math.h>

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
    double scan_pos_x;     /* the scan position in x */
    double scan_pos_z;     /* the scan position in z */
    int make_sphere;       /* Should the analytic sphere be added to the model */
    double *dsphere_c;     /* The pinhole-plate sample distance (for use with 
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
    int32_t *cntr_detected;       /* The number of detected rays */
    int killed;              /* The number of killed rays */
    int32_t *numScattersRay; /* The number of sample scatters that each
                              * ray has undergone */
    
    /* Declare other variables */
    int i;
    gsl_rng *my_rng;
    int sample_index, sphere_index, plate_index;
    int detector;
    
    /* Declare structs */
    Surface3D Sample;
    NBackWall Plate;
    AnalytSphere the_sphere;
    
    /*******************************************************************************/
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 21) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Twenty one inputs required for tracingMex.");
    }
    if (nlhs != 3) {
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
    V = mxGetPr(prhs[0]);
    F = mxGetPr(prhs[1]); /* Reading in as double not int - cast later in code */
    N = mxGetPr(prhs[2]);
    C = mxGetPr(prhs[3]);
    P = mxGetPr(prhs[4]);
    maxScatters = (int)mxGetScalar(prhs[5]); /* mxGetScalar gives a double */
    scan_pos_x = mxGetScalar(prhs[6]);
    scan_pos_z = mxGetScalar(prhs[7]);
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
    nrays = (int)mxGetScalar(prhs[18]);
    source_model = (int)mxGetScalar(prhs[19]);
    source_parameters = mxGetPr(prhs[20]);
    ntriag_sample = mxGetN(prhs[1]);
    
    
    /**************************************************************************/
    
    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;
    
    /* Set up the GSL random number generator */
    my_rng = setupGSL();
    
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
    the_sphere = set_up_sphere(make_sphere, scan_pos_x, scan_pos_z, sphere_c,
        sphere_r, sphere_diffuse, sphere_parameters, sphere_index);
        
    /* 
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the 
     * difference in indexing between MATLAB and C.
     */
    plhs[0] = mxCreateNumericMatrix(1, n_detector, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, n_detector*maxScatters, mxINT32_CLASS, mxREAL);
    
    /* Pointers to the output matrices so we may change them*/
    cntr_detected = (int32_t*)mxGetData(plhs[0]);
    numScattersRay = (int32_t*)mxGetData(plhs[2]);
    
    /**************************************************************************/
    
    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (i = 0; i < nrays; i++) {
        Ray3D the_ray;
        int detected;
        
        the_ray = create_ray_source(source_parameters[0], &source_parameters[1], 
            source_parameters[3], source_parameters[4], source_model, my_rng, 
            source_parameters[5]);
        
        detected = trace_ray_simpleMulti(&the_ray, &killed, cntr_detected,
            maxScatters, Sample, Plate, the_sphere, my_rng, &detector);
                
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
    
    /* Free the space used by the random number generator */
    gsl_rng_free(my_rng);
    
    return;
}


