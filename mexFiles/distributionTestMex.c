/*
 * Copyright (c) 2020-21, Dan Seremet, Sam Lambrick.
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
#include "extract_inputs.h"
#include "atom_ray_tracing3D.h"

void print_material_mex(Material const * const mat);

/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    const int N_INPUTS = 5;
    const int N_OUTPUTS = 3;
    
    /* Input variables */
    int n_rays;
    double *direction;
    double *normal;
    double * lattice;
    Material material;
    
    /* Output variables */
    double * thetas;
    double * phis;
    double * final_dir;
    
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
    lattice = mxGetDoubles(prhs[4]);
    
    // check and extract material properties
    char *name[] = {"default"};
    get_materials(prhs[2], name, &material);

    // print parameters
    mexPrintf("Running distribution_test_mex with:\n");
    mexPrintf("n_rays = %d \t", n_rays);
    mexPrintf("direction = [%.2f %.2f %.2f] \t", direction[0], direction[1],
              direction[2]);
    mexPrintf("normal = [%.2f %.2f %.2f] \n", normal[0], normal[1], normal[2]);
    print_material_mex(&material);
    mexPrintf("\n");

    // DONE extracting params. Setup RNG and proceed to calculation
    mexPrintf("distribution_test.c -- inputs successfully extracted\n");
    
    /* Seed the random number generator with the current time */
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
    
    /* Set up the MTwister random number generator */
    seedRand(t, &myrng);

    // allocate output array and get pointer
    plhs[0] = mxCreateDoubleMatrix(n_rays, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n_rays, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n_rays, 3, mxREAL);
    thetas = mxGetDoubles(plhs[0]);
    phis = mxGetDoubles(plhs[1]);
    final_dir = mxGetDoubles(plhs[2]);

    /* project the initial direction onto the surface plane to find azimuth angle phi
     * this looks like d' = d - (n.d)n */
    double tmp;
    dot(direction, normal, &tmp);
    propagate(direction, normal, -tmp, dir_projection);
    if ((dir_projection[0] + dir_projection[1] + dir_projection[2]) < 1e-6) {
        dir_projection[0] = 1;
        dir_projection[1] = 1;
        dir_projection[2] = 1;
    }
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

        material.func(normal, lattice, direction, new_dir, material.params, &myrng);
        normalise(new_dir);
        double tmp;
        dot(new_dir, normal, &tmp);
        propagate(new_dir, normal, -tmp, new_dir_proj);
        normalise(new_dir_proj);

        dot(new_dir, normal, &tmp);
        thetas[i] = acos(tmp);

        dot(new_dir_proj, perpendicular, &sin_phi);
        dot(new_dir_proj, dir_projection, &cos_phi);
        phis[i] = atan2(sin_phi, cos_phi);
        
        for (int j = 0; j < 3; j++)
            final_dir[3*i+j] = new_dir[j];
    }

    mexPrintf("done.\n\n");
    return;
}

/* print details of Material struct */
void print_material_mex(Material const * const mat) {
    mexPrintf("\tMAT %-10s func %-12s", mat->name, mat->func_name);
    for (int i = 0; i < mat->n_params; i++)
        mexPrintf(" %.2f ", mat->params[i]);
}
