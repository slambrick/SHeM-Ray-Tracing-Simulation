/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Tests the sampling fo the cosine distribution used in the SHeM simulation. 
 * Could be used as a template to test other distributions.
 *
 * The calling syntax is:
 *
 *  positions = test_cosineScatter(N)
 *  
 *  INPUTS:
 *   N - int, the number of directions to generate
 * 
 *  OUTPUTS:
 *   positions - A 2d, double, array of directions generated accoring to a 
 *               normally centred cosine distribution.
 * 
 * This is a MEX file for MATLAB.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mex.h>
#include <gsl/gsl_rng.h>
#include "small_functions.h"

/* 
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    
    /**************************************************************************/
    
    /* Declare the input variables */
    int N;
    double* normal;
    
    /* Declare the output variables */
    double *some_positions;
    
    /* Declare other variables */
    int i;
    gsl_rng *myrng;
    double init_dir[3];
    
    init_dir[0] = 1/sqrt(2);
    init_dir[1] = -1/sqrt(2);
    init_dir[2] = 0;
    
    /**************************************************************************/
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:test_cosineScatter:nrhs", 
                          "One inputs required for test_cosineScatter.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:test_cosineScatter:nrhs", 
                          "One output required for test_cosineScatter.");
    }
    
    /* Checks the correct type of the input. */
    if (!mxIsDouble(prhs[0]) || 
        mxIsComplex(prhs[0]) || 
        mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:test_cosineScatter:notScalar",
                          "Input must be a real scalar.");
    }
    
    /* Checks that the number provided is positive */
    if (prhs[0] <= 0) {
        mexErrMsgIdAndTxt("MyToolbox:test_cosineScatter:badInput",
                          "Must specify a positive number of samples to take");
    }
    
    /**************************************************************************/
    
    /* Read the input variables */
    N = (int)mxGetScalar(prhs[0]);
    normal = mxGetPr(prhs[1]);
    
    /* Pointers to output variables */
    plhs[0] = mxCreateDoubleMatrix(3, N, mxREAL);
    some_positions =  mxGetPr(plhs[0]);
    
    /* Seed the random number generator. */
    myrng = setupGSL();
    
    for (i = 0; i < N; i++) {
        int k;
        double random_dir[3];
        
        cosineSpecularScatter(normal, init_dir, random_dir, myrng);
        
        for (k = 0; k < 3; k++) {
            int n;
            n = 3*i + k;
            some_positions[n] = random_dir[k];
        }
    }
    
    return;
}


