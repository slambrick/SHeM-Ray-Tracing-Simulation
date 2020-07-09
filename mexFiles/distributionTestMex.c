/*
 * Copyright (c) 2020, Dan Seremet, Sam Lambrick.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * MEX function for testing what the probability distributions do.
 * It produced the angle between the reflected direction and the normal.
 *
 * The calling syntax is:
 * angles = distribution_test(n_rays, direction, material, normal);
 *
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
#include "small_functions3D.h"
#include "ray_tracing_structs3D.h"

/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    const int N_INPUTS = 4;
    const int N_OUTPUTS = 2;
    
    /* Input variables */
    int n_rays;
    double *direction;
    double *normal;
    Material material;
    
    /* Output variables */
    double * thetas;
    double * phis;
    
    /* Other varibles */
    double dir_projection[3];
    double perpendicular[3];
    int i;
    
    /* For random number generation */
    struct timeval tv;
    unsigned long t;
    MTRand myrng;
    
    /* Check for the right number of inputs and outputs */
    if (nrhs != N_INPUTS)
        mexErrMsgIdAndTxt("test:distribution:nrhs",
                          "%d inputs required", N_INPUTS);
    if (nlhs != N_OUTPUTS)
        mexErrMsgIdAndTxt("test:distribution:nrhs",
                          "%d outpus required", N_OUTPUTS);
    
    /* Check for the right type of inputs */
    if(!mxIsScalar(prhs[0]))
        mexErrMsgIdAndTxt("test:distribution:n_rays",
                          "n_rays should be a scalar");

    if(!mxIsDouble(prhs[1]) || mxGetN(prhs[1]) != 3)
        mexErrMsgIdAndTxt("test:distribution:direction",
                          "initial ray direction should be 3-vector");
    
    if(!mxIsDouble(prhs[3]) || mxGetN(prhs[3]) != 3)
        mexErrMsgIdAndTxt("test:distribution:normal",
                          "initial ray normal should be 3-vector");
    
    /* Get the input variables */
    n_rays = (int)mxGetScalar(prhs[0]);
    direction = mxGetDoubles(prhs[1]);
    normal = mxGetDoubles(prhs[3]);
    
    // check and extract material properties
    char *name[] = {"default"};
    get_materials(prhs[2], name, &material);

    // print parameters
    mexPrintf("Running distribution_test_mex with:\n");
    mexPrintf("n_rays = %d \t", n_rays);
    mexPrintf("direction = [%.2f %.2f %.2f] \t", direction[0], direction[1],
              direction[2]);
    mexPrintf("normal = [%.2f %.2f %.2f] \n", normal[0], normal[1], normal[2]);
    print_material(&material);
    mexPrintf("\n");

    // DONE extracting params. Setup RNG and proceed to calculation
    mexPrintf("distribution_test.c -- inputs successfully extracted\n");
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    /* Set up the MTwister random number generator */
    myrng = seedRand(t);

    // allocate output array and get pointer
    plhs[0] = mxCreateDoubleMatrix(n_rays, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_rays, 1, mxREAL);
    thetas = mxGetDoubles(plhs[0]);
    phis = mxGetDoubles(plhs[1]);

    /* project the initial direction onto the surface plane to find azimuth angle phi
     * this looks like d' = d - (n.d)n */
    propagate(direction, normal, -dot(direction, normal), dir_projection);
    normalise(dir_projection);

    /* also calculate perpendicular -- needed to find sign of phi
     * this is p = n x d' */
    cross(normal, dir_projection, perpendicular);

    // another way to do this is (uncomment):
    // cross(normal, direction, perpendicular);
    // cross(perpendicular, normal, dir_projection);

    mexPrintf("Scattering your rays... ");
    // for the specified number of rays, perform the scattering event
    // and accumulate the angle of (new_dir, normal) in the output array
    
    for (i = 0; i < n_rays; i++) {
        double cos_phi, sin_phi;
        double new_dir_proj[3];
        double new_dir[3] = {0, 1, 0};

        material.func(normal, direction, new_dir, material.params, &myrng);
        normalise(new_dir);
        propagate(new_dir, normal, -dot(new_dir, normal), new_dir_proj);
        normalise(new_dir_proj);

        thetas[i] = acos(dot(new_dir, normal));

        sin_phi = dot(new_dir_proj, perpendicular);
        cos_phi = dot(new_dir_proj, dir_projection);
        phis[i] = atan2(sin_phi, cos_phi);
    }

    mexPrintf("done.\n\n");
    return;
}


