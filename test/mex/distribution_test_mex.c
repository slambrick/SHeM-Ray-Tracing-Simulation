/*
 * Copyright (c) 2020, Dan Seremet.
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

#include <gsl/gsl_rng.h>
#include <math.h>

#include "extract_inputs.h"
#include "common_helpers.h"
#include "small_functions3D.h"
#include "ray_tracing_structs3D.h"
#include "distributions.h"


/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    const int N_INPUTS = 4;
    const int N_OUTPUTS = 2;

    /* Check for the right number of inputs and outputs */
    if (nrhs != N_INPUTS)
        mexErrMsgIdAndTxt("test:distribution:nrhs",
                          "%d inputs required", N_INPUTS);
    if (nlhs != N_OUTPUTS)
        mexErrMsgIdAndTxt("test:distribution:nrhs",
                          "%d outpus required", N_OUTPUTS);

    // check and extract number of rays
    if(!mxIsScalar(prhs[0]))
        mexErrMsgIdAndTxt("test:distribution:n_rays",
                          "n_rays should be a scalar");
    int n_rays = (int)mxGetScalar(prhs[0]);

    // check and extract initial direction of rays
    if(!mxIsDouble(prhs[1]) || mxGetN(prhs[1]) != 3)
        mexErrMsgIdAndTxt("test:distribution:direction",
                          "initial ray direction should be 3-vector");
    double * direction = mxGetDoubles(prhs[1]);

    // check and extract material properties
    Material material; char *name[] = {"default"};
    get_materials(prhs[2], name, &material);

    // check and extract surface normal
    if(!mxIsDouble(prhs[3]) || mxGetN(prhs[3]) != 3)
        mexErrMsgIdAndTxt("test:distribution:normal",
                          "initial ray normal should be 3-vector");
    double * normal = mxGetDoubles(prhs[3]);

    // print parameters
    mexPrintf("Running distribution_test_mex with:\n");
    mexPrintf("n_rays = %d \t", n_rays);
    mexPrintf("direction = [%.2f %.2f %.2f] \t", direction[0], direction[1], direction[2]);
    mexPrintf("normal = [%.2f %.2f %.2f] \n", normal[0], normal[1], normal[2]);
    print_material(&material);
    mexPrintf("\n");

    // DONE extracting params. Setup RNG and proceed to calculation
    mexPrintf("distribution_test.c -- inputs successfully extracted\n");
    gsl_rng * my_rng = setupGSL();

    // allocate output array and get pointer
    plhs[0] = mxCreateDoubleMatrix(n_rays, 1, mxREAL);
    double * thetas = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(n_rays, 1, mxREAL);
    double * phis = mxGetDoubles(plhs[1]);

    // project the initial direction onto the surface plane to find azimuth angle phi
    // this looks like d' = d - (n.d)n
    double dir_projection[3], perpendicular[3];
    propagate(direction, normal, -dot(direction, normal), dir_projection);
    normalise(dir_projection);

    // also calculate perpendicular -- needed to find sign of phi
    // this is p = n x d'
    cross(normal, dir_projection, perpendicular);

    // another way to do this is (uncomment):
    // cross(normal, direction, perpendicular);
    // cross(perpendicular, normal, dir_projection);

    mexPrintf("Reflecting your rays...\n");
    // for the specified number of rays, perform the scattering event
    // and accumulate the angle of (new_dir, normal) in the output array
    for(int iray = 0; iray < n_rays; iray++) {
        double new_dir[3] = {0, 1, 0};

        material.func(normal, direction, new_dir, material.params, my_rng);
        normalise(new_dir);

        double new_dir_proj[3];
        propagate(new_dir, normal, -dot(new_dir, normal), new_dir_proj);
        normalise(new_dir_proj);

        thetas[iray] = acos(dot(new_dir, normal));

        double sin_phi = dot(new_dir_proj, perpendicular);
        double cos_phi = dot(new_dir_proj, dir_projection);
        phis[iray] = atan2(sin_phi, cos_phi);
    }

    mexPrintf("\n");
    return;
}


