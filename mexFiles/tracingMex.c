/*
 * Copyright (c) 2018, Sam Lambrick.
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
 *                 dist_to_sample, sphere_r, diffuse);
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
#include <gsl/gsl_rng.h>
#include "small_functions.h"
#include "tracing_functions.h"
#include <math.h>

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
    int nrays;             /* number of rays */
    int ntriag_sample;     /* number of sample triangles */
    int maxScatters;       /* Maximum number of scattering events per ray */
    double *VS;            /* pinhole plate triangle vertices */
    double *FS;            /* pinhole plate triangle faces */
    double *NS;            /* pinhole plate triangle normals */
    double *CS;            /* pinhole plate triangle diffuse level*/
    int ntriag_plate;      /* number of pinhole plate triangles */
    double *backWall;      /* y-coordinate of the back of the pinhole plate, the 
                            * depth in x and then z of the pinhole plate */
    double scan_pos_x;     /* the scan position in x */
    double scan_pos_z;     /* the scan position in z */
    int make_sphere;       /* Should the analytic sphere be added to the model */
    double dist_to_sample; /* The pinhole-plate sample distance (for use with 
                            * the analytic sphere) */
    double sphere_r;       /* Radius of the analytic sphere if it to be made */
    double sphere_diffuse; /* The scattering off of the analytic sphere */
    
    /* Declare the output variables */
    int cntr;               /* The number of detected rays */
    int killed;             /* The number of killed rays */
    double *final_pos;      /* The final positions of the detected rays */
    double *final_dir;      /* The final directions of the detected rays */
    double *numScattersRay; /* The number of sample scatters that each detected
                             * ray has undergone */
    
    /* Declare other variables */
    int i;
    gsl_rng *myrng;
    
    /*******************************************************************************/
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != 18) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Eighteen inputs required for tracingMex.");
    }
    if (nlhs != 5) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs", 
                          "Five outpus required for tracingMex.");
    }
    
    /* Check the type of the inputs */
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
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
    }
    
    /**************************************************************************/
    
    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    ray_pos = mxGetPr(prhs[0]);
    ray_dir = mxGetPr(prhs[1]);
    V = mxGetPr(prhs[2]);
    F = mxGetPr(prhs[3]); /* Reading in as double not int - cast to int later in code */
    N = mxGetPr(prhs[4]);
    C = mxGetPr(prhs[5]);
    maxScatters = (int)mxGetScalar(prhs[6]); /* mxGetScalar gives a double */
    VS = mxGetPr(prhs[7]);
    FS = mxGetPr(prhs[8]);
    NS = mxGetPr(prhs[9]);
    CS = mxGetPr(prhs[10]);
    backWall = mxGetPr(prhs[11]);
    scan_pos_x = mxGetScalar(prhs[12]);
    scan_pos_z = mxGetScalar(prhs[13]);
    make_sphere = (int)mxGetScalar(prhs[14]); /* mxGetScalar gives a double */
    dist_to_sample = mxGetScalar(prhs[15]);
    sphere_r = mxGetScalar(prhs[16]);
    sphere_diffuse = mxGetScalar(prhs[17]);
    nrays = mxGetN(prhs[0]);
    ntriag_sample = mxGetN(prhs[3]);
    ntriag_plate = mxGetN(prhs[8]);
    
    /**************************************************************************/
    
    /* 
     * Number of rays that enter the detector and those that go into the detector
     * after a single scatter 
     */
    cntr = 0;
    
    /* Number of rays that are killed as they have scattered too many times */
    killed = 0;
    
    /* 
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the 
     * difference in indexing between MATLAB and C.
     */
    plhs[2] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(3, nrays, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, nrays, mxREAL);

    
    /* Pointers to the ouput matrices so we may change them*/
    final_pos       = mxGetPr(plhs[2]);
    final_dir       = mxGetPr(plhs[3]);
    numScattersRay  = mxGetPr(plhs[4]);
    
    /* Set up the GSL random number generator */
    myrng = setupGSL();
    
    /**************************************************************************/
    
    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    for (i = 0; i < nrays; i++) {
        int dead;
        double old[3];	/* The previous location the ray was at */
        double e[3];
        double d[3];
        int k, nScatters;
        int current_tri;
        int current_surface;
        int n_allScatters;
        
        /* 
         * The total number of scattering events undergone (sample and pinhole 
         * plate) 1000 events are allowed in total. A seperate limit is placed
         * on the number of scattering events off of the sample.
         */
        n_allScatters = 0;
        
        /* A method for storing the current triangle the ray is on, -1 is not 
         * on any 
         */
        current_tri = -1;
        
        /* A method for storing which surface the ray is on 0=sample, 1=plate, 
         * -1=none, -2=sphere 
         */
        current_surface = -1;
        
        /*
         * If the ray is 'dead' it is no longer of interest. It may be dead in two 
         * ways, either by not meeting any surface of by going into the detector. 
         *      dead = 1, it did not meet any surfaces 
         *      dead = 2, it went into the detector
         */
        dead = 0;
        
        /* A counter for the number of scattering events off of the sample */
        nScatters = 0;
        
        /* The position and direction of the ray
         * The MATLAB matrices are imported as 1D arrays in C, hence the
         * indexing 
         */
        for (k = 0; k < 3; k++) {
            int n;
            n = 3*i + k;
            e[k] = ray_pos[n];
            d[k] = ray_dir[n];
        }
        
        /* 
         * Keep propogating the ray until it is deemed 'dead', by either not 
         * intersecting either the sample or the pinhole plate. We count the 
         * number of scattering events all rays undergo storing in 
         * numScattersRay, we add 100 extra if it is 'killed' for exceeding the
         * allowed number of scatters.
         */
        while (!dead) {
            double new_dir[3];
            double middle[3];
            int yes;
            
            /* What is yes? */
            yes = 0;
            
            /* The ray is dead unless we hit something */
            dead = 1;
            
            if ((nScatters > maxScatters) || (n_allScatters > 1000)) {
                /* Ray has exceeded the maximum number of scatters, kill it */
                numScattersRay[i] = nScatters + 100;
                killed += 1;
                break;
            }
            
            /******************************************************************/
            /* 
             * Try to scatter of sample. This only tries to scatter off of the 
             * sample and not the pinhole plate, must be done first as the newly
             * generated rays will in about half of the cases hit the inside of
             * the pinhole plate, which will cause weird (and very wrong) 
             * results. Needs only be done when the ray has just been generated.
             * 
             * dead = 1 by scatterOffSurface -- did not hit the sample
             * dead = 0 by scatterOffSurface -- scattered off the sample
             * 
             * If the ray has not hit the sample then it is immediatlly dead.
             */
            if (nScatters == 0) {
                dead = scatterOffSurface(e, d, ntriag_sample, V, N, F, C, myrng, 
                    &current_tri, &current_surface, scan_pos_x, scan_pos_z, 
                    make_sphere, dist_to_sample, sphere_r, sphere_diffuse);
                
                if (dead == 0) {
                    /* Hit the sample */
                    nScatters++;
                    n_allScatters++;
                } else {
                    /* Did not hit the sample, therefore the ray is dead */
                    numScattersRay[i] = nScatters;
                    
                    /* Move onto the next ray */
                    continue;
                }
            }
            
            /* Try to scatter of both surfaces. All but the first iteration */
            dead = scatterSurfaces(e, d, ntriag_sample, ntriag_plate, V, N, F, 
                C, VS, NS, FS, CS, myrng, &current_tri, &current_surface, 
                backWall, scan_pos_x, scan_pos_z, make_sphere, dist_to_sample, 
                sphere_r, sphere_diffuse);
            
            /******************************************************************/
            /* Update counters */
            
            switch (dead) {
                case 2:
                    /* Detected */
                    for (k = 0; k < 3; k++) {
                        int n;
                        n = 3*i + k;
                        final_pos[n] = e[k];
                        final_dir[n] = d[k];
                    }
                    cntr += 1;
                    numScattersRay[i] = nScatters;
                    break;
                case 1:
                    /* Did not hit a surface or get detected */
                    numScattersRay[i] = nScatters;
                    break;
                case 0:
                    /* Hit a surface */
                    n_allScatters++;
                    
                    /* Hit the sample */
                    if ((current_surface == 0) || current_surface == -2) {
                        nScatters++;
                    }
                    break;
            }
        }
    }
    
    /* Output number of rays went into the detector */
    plhs[0] = mxCreateDoubleScalar(cntr);
    plhs[1] = mxCreateDoubleScalar(killed);
    
    /* Free the space used by the random number generator */
    gsl_rng_free(myrng);
    
    return;
}


