/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the  
 * GNU/GPL-3.0-or-later.
 */
#include "mex.h"
#include "intersect_detection2D.h"
#include <sys/time.h>
#include <stdlib.h>
#include "mtwister.h"
#include <stdint-gcc.h>

/******************************************************************************/
/*                The gateway function for calling from MATLAB                */
/******************************************************************************/

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    
    /**************************************************************************/
    /* Declare variables */
    
    /* Declare input variables */
    double *vertices;
    double *normals;
    double *ray_pos;
    double *ray_dir;
    int nrays;
    int nelements;
    int scattering;
    double *scattering_parameters;
    
    /* Declare output variables */
    double *final_dir;
    double *final_pos;
    mxInt32 *num_scatters;
    
    /* Declare other variables */
    Surface2D Sample;
    Ray2D *rays;
    Rays2D all_rays;
    int i;
    struct timeval tv;
    unsigned long t;
    MTRand my_rng;
    
    /**************************************************************************/
    /* Checks on inputs */
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MyToolbox:scatterRaysMex:nrhs", 
                          "Six inputs required for tracingMex.");
    }
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:scatterRaysMex:nlhs", 
                          "Three outputs required for tracingMex.");
    }
    
    /**************************************************************************/
    /* Set up memory, pointers to input and output variables, allocation of 
     * memory etc. */
    
    /* Read the input variables */
    vertices = mxGetPr(prhs[0]);
    normals = mxGetPr(prhs[1]);
    nelements = mxGetN(prhs[1]);
    ray_pos = mxGetPr(prhs[2]);
    nrays = mxGetN(prhs[2]);
    ray_dir = mxGetPr(prhs[3]);
    scattering = (int)mxGetScalar(prhs[4]);
    scattering_parameters = mxGetPr(prhs[5]);
    
    /* Create the output matrices and assign a pointer to it */
    plhs[0] = mxCreateDoubleMatrix(2, nrays, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(2, nrays, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, nrays, mxINT32_CLASS, mxREAL);
    final_dir = mxGetPr(plhs[0]);
    final_pos = mxGetPr(plhs[1]);
    num_scatters = (mxInt32*)mxGetData(plhs[2]);
    
    /* TODO: Move to Dan's way of doing scattering. */
    
    /* Put the surface into a struct */
    Sample.V = vertices;
    Sample.N = normals;
    Sample.n_elements = nelements;
    Sample.scattering = scattering;
    Sample.scattering_parameters = scattering_parameters;
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    srand(t);
    /* Set up the MTwister random number generator */
    my_rng = seedRand(t);
    Sample.my_rng = &my_rng;
    
    /* Allocate appropriate memory to the array of Ray structs */
    rays = malloc(nrays * sizeof(*rays));
    
    /* Put the rays into an array of structs */
    compse_rays2D(rays, ray_pos, ray_dir, nrays);
    
    /* Put the rays and the number of rays into a single Rays struct */
    all_rays.rays = rays;
    all_rays.nrays = nrays;
    
    /**************************************************************************/
    /* Perform the simulation */
    
    /* Call the main function that performs the computation */
    scatterRays2D(Sample, all_rays);
    
    /**************************************************************************/
    /* Put the results in output variables and deallocate memory */
    
    /* Put the final directions into the output variable */
    get_directions2D(all_rays, final_dir);
    
    /* Put the final positions into the output variable */
    get_positions2D(all_rays, final_pos);
    
    /* Put the number of scattering events in the output variabel */
    for (i = 0; i < all_rays.nrays; i++) {
        num_scatters[i] = all_rays.rays[i].nscatters;
    }
    
    /* De-allocate the memory assigned by malloc and the random number generator */
    
    free(rays);
    
    return;
}




