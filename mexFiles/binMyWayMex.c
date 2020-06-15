/*
 * Copyright (c) 2018, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the 
 * GNU/GPL-3.0-or-later.
 *
 * 
 * Takes a 1D vector of positive integers and counnts the number of occurances
 * of each number, counts for all values from 1 upto the maximum numbe
 * which is given. Any values above the maximum number are ignored.
 * 
 * The calling syntax is:
 *  
 *  histRay = binMyWayMex(numScattersRay, maxScatters)
 * 
 *  INPUTS:
 *   numScattersRay - 1D positive int column array to bin, must have more than
 *                    two elements
 *   maxScatters    - The number to bin up to, must be an integer greater than 1
 * 
 *  OUTPUTS
 *   histRay - 1D row array of the binned variables
 * 
 *  EXAMPLE:
 *   % The call
 *   binMyWayMex([2,5,6,2,4,6,5,5,2], 5);
 *   % would give the output
 *   [0, 3, 0, 1, 3].
 * 
 * This is a MEX file for MATLAB.
 */
#include "mex.h"

void binIt(int numScattersRay[], int histRay[], int maxScatters, int nRays);

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {
    
    /* Declare input variables */
    int *numScattersRay;
    int maxScatters;
    int nRays;
    
    /* Declare output variables */
    int *histRay;
    
    /* Declare other variables */
    int i;
    
    /*************************************************************************/
    
    /* Check for the right number of inputs */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:nrhs", 
                          "Two inputs required for binMyWayMex.");
    }
    
    /* Check for the right number of outputs */
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:nlhs", 
                          "One output required for binMyWayMex.");
    }
    
    /* Check for the right input types*/
    if (!mxIsDouble(prhs[1]) || 
        mxIsComplex(prhs[1]) || 
        mxGetNumberOfElements(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:notScalar", 
                          "Number of rays must be a real scaler.");
    }
    if (!mxIsInt32(prhs[0]) || mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:notDouble", 
                          "Input vector must be of type int.");
    }
    
    /* More checking on the inputs */
    if (mxGetScalar(prhs[1]) < 1) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:badInput",
                          "The maximum bin value must be greater than 1.");
    }
    if (mxGetN(prhs[0]) <= 2) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:badInput",
                          "Input array needs to have more than two rows.");
    }
    if (mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:badInput",
                          "Input array should be a column vector."); 
    }
    
    /*************************************************************************/
    
    /* Read input variables */
    numScattersRay = (int*)mxGetData(prhs[0]);
    maxScatters = (int)mxGetScalar(prhs[1]);
    nRays = mxGetN(prhs[0]);
    
    /* Create output varibles */
    plhs[0] = mxCreateNumericMatrix(1, maxScatters, mxINT32_CLASS, mxREAL);
    histRay = (int*)mxGetData(plhs[0]); /* Safe to cast */
    
    /* Check that non of the numbers in the input array are negative */
    for (i = 0; i < nRays; i++) {
        if (numScattersRay[i] < 0) {
            mexErrMsgIdAndTxt("MyToolbox:binMyWayMex:badInput",
                "Input array has negative numbers.");
        }
    }
    
    /*************************************************************************/
    
    /* Numerical routine */
    binIt(numScattersRay, histRay, maxScatters, nRays);
    
    return;
}

void binIt(int numScattersRay[], int histRay[], int maxScatters, int nRays) {
    int i, j;
    
    for (i = 0; i < nRays; i++) {
        for (j = 1; j <= maxScatters; j++) {
            if (j == numScattersRay[i]) {
                histRay[j - 1] += 1;
            }
        }
    }
}

