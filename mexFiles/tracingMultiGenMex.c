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

#include <mex.h>
#include <matrix.h>

#include <gsl/gsl_rng.h>
#include <stdint.h>
#include <math.h>

// #include "trace_ray.h"
// #include "small_functions3D.h"
#include "common_helpers.h"
#include "ray_tracing_structs3D.h"

/*
 * Take the elements from a MATLAB cell array of strings
 * and put them in a C matrix of characters (or array of pointers to strings, really).
 * Returns how many strings were extracted.
 *
 * NB: this function will not allocate space for the character matrix.
 * You should allocate space in the same function that frees it, i.e. the caller.
 */
int getStrCellArray(const mxArray * cell_array, char ** strings) {
    int num = mxGetN(cell_array);
    mxArray * cell;

    for(int icell = 0; icell < num; icell++) {
        cell = mxGetCell(cell_array, icell);
        if (!mxIsChar(cell))
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:strings",
                              "Each cell must be a char array.");
        if (mxGetM(cell) != 1)
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:strings",
                              "Each cell must be a row vector.");

        strings[icell] = mxArrayToString(cell);
    }
    return num;
}

/*
 * Extract the sphere C struct from the MATLAB struct array.
 * INPUTS: theSphere = mxArray containing ONE sphere struct
 *         index = the surface_index of the sphere surface
 */
AnalytSphere getSphere(const mxArray * theSphere, int index) {
    // check if sphere is one struct
    if(!mxIsStruct(theSphere))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:theSphere",
                          "Must be struct array");
    if(mxGetN(theSphere) != 1 || mxGetM(theSphere) != 1)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:theSphere",
                          "Must be one struct");

    // get the centre
    mxArray * field = mxGetField(theSphere, 0, "c");
    if (!mxIsDouble(field) || mxGetN(field) != 3)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:theSphere",
                          "Centre must be array of 3 doubles.");
    double * sphere_c = mxGetDoubles(field);

    // get the radius - one double
    field = mxGetField(theSphere, 0, "r");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:theSphere",
                          "Sphere radius must be scalar");
    double sphere_r = mxGetScalar(field);

    // get the 'make' - one integer
    field = mxGetField(theSphere, 0, "make");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:theSphere",
                          "Sphere make must be scalar");
    int make_sphere = (int)mxGetScalar(field);

    // now extract the material -- a nested struct
    mxArray * sph_material = mxGetField(theSphere, 0, "material");
    if(!mxIsStruct(sph_material))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:sphereMaterial",
                          "Must be struct array");
    if(mxGetN(sph_material) != 1 || mxGetM(sph_material) != 1)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:sphereMaterial",
                          "Must be one struct");
    int num_params = mxGetN(mxGetField(sph_material, 0, "params"));
    double * params = mxGetDoubles(mxGetField(sph_material, 0, "params"));
    char * function = mxArrayToString(mxGetField(sph_material, 0, "function"));

    Material mat = set_up_material("sphere", function, params, num_params);

    AnalytSphere sph = set_up_sphere(make_sphere, sphere_c, sphere_r, mat, index);
    return sph;
}

/*
 * Print faces of the sample surface to the console, as received from the MATLAB mex func.
 */
void printSample(const int num_faces, const int32_t * F, const double * N, char ** C) {
    for(int iface = 0; iface < num_faces; iface++) {
        mexPrintf("\n\t FACE %2d", iface);
        mexPrintf("\t V %3d %3d %3d", F[3*iface], F[3*iface+1], F[3*iface+2]);
        mexPrintf("\t N %+f %+f %+f", N[3*iface], N[3*iface+1], N[3*iface+2]);
        mexPrintf("\t C %s", C[iface]);
    }
    mexPrintf("\n");
}


/**
 * Package the names, function names and parameters in the input into
 * an array of Material pointers. Return the number of materials found.
 */
int getMaterials(const mxArray * names, const mxArray * functions,
                 const mxArray * params, Material * materials) {

    if(!mxIsCell(names) || !mxIsCell(functions) || !mxIsCell(params))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:materials",
                          "Materials arrays must be cell arrays.");

    int num = mxGetN(names);
    if(mxGetN(functions) != num || mxGetN(params) != num)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:materials",
                          "Dimensions of materials arrays must be the same.");

    // use cell array extraction for names and functions
    char ** names_extr = mxCalloc(num, sizeof(char*));
    char ** functions_extr = mxCalloc(num, sizeof(char*));

    getStrCellArray(names, names_extr);
    getStrCellArray(functions, functions_extr);

    for(int idx = 0; idx < num; idx++) {
        mxArray * params_cell = mxGetCell(params, idx);
        int num_params = mxGetN(params_cell);

        mexPrintf("\nNum params %d", num_params);

        if(!mxIsDouble(params_cell))
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:materials",
                              "Parameters of each material must be double array.");
        double * params_extr = mxGetDoubles(params_cell);
        materials[idx] = set_up_material(
            names_extr[idx], functions_extr[idx], params_extr, num_params);
    }

    mxFree(names_extr); mxFree(functions_extr);
    return num;
}


/*
 * The gateway function.
 * lhs = left-hand-side, outputs
 * rhs = right-hand-side, inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {


    const int NINPUTS = 13;
    const int NOUTPUTS = 3;

    /* Declare the input variables */
    int ntriag_sample;     /* number of sample triangles */
    double *V;             /* sample triangle vertices 3xn */
    int32_t *F;            /* sample triangle faces 3xM */
    double *N;             /* sample triangle normals 3xM */
    char **C;              /* sample material keys, length M */
    Material *M;           /* materials of the sample */
    int n_rays;            /* number of rays */
    int maxScatters;       /* Maximum number of scattering events per ray */
    int source_model;
    double *source_parameters;

    mexPrintf("\n\n Now in tracingMultiGenMex \n\n");

    /* Declare the output variables */
    int32_t * cntr_detected;       /* The number of detected rays */
    int killed = 0;              /* The number of killed rays */
    int32_t * numScattersRay; /* The number of sample scatters that each
                               * ray has undergone */

    /* Declare other variables */
    int detector;
    int status;         // a status indicator for functions returning success/error codes
    mxArray * field;    // an array pointer for unpacking nested cell arrays or structs

    Surface3D sample;
    NBackWall plate;
    AnalytSphere sphere;

    /* Indexing the surfaces, -1 refers to no surface */
    int sample_index = 0, plate_index = 1, sphere_index = 2;

    /* Check for the right number of inputs and outputs */
    if (nrhs != NINPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "%d inputs required for tracingMultiGenMex.", NINPUTS);
    }
    if (nlhs != NOUTPUTS) {
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:nrhs",
                          "%d outpus required for tracingMultiGenMex.", NOUTPUTS);
    }

    /* Read the input variables.
     * NOTE: mxGetScalar always returns type double. In cases that the input in
     *       MATLAB were of type int it is safe to cast from double to int here.
     */
    V = mxGetDoubles(prhs[0]);

    ntriag_sample = mxGetN(prhs[1]);
    F = mxGetInt32s(prhs[1]);
    N = mxGetDoubles(prhs[2]);

    // read in the material keys
    C = mxCalloc(ntriag_sample, sizeof(char*));
    getStrCellArray(prhs[3], C);
    // printSample(ntriag_sample, F, N, C);

    // get the sphere from struct
    sphere = getSphere(prhs[4], sphere_index);

    // extract plate properties from thePlate cell array containing plate options
    const mxArray * plate_opts = prhs[5];
    if(!mxIsCell(plate_opts))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate_opts",
                          "Must be cell array");
    if(mxGetN(plate_opts) != 5 && mxGetM(plate_opts) != 5)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate_opts",
                          "Must be of length 5");
    plate.plate_represent = (int)mxGetScalar(mxGetCell(plate_opts, 0));
    plate.n_detect = (int)mxGetScalar(mxGetCell(plate_opts, 1));
    plate.circle_plate_r = mxGetScalar(mxGetCell(plate_opts, 2));
    plate.aperture_axes = mxGetDoubles(mxGetCell(plate_opts, 3));
    plate.aperture_c = mxGetDoubles(mxGetCell(plate_opts, 4));
    plate.surf_index = plate_index;

    // materials
    int num_materials = mxGetN(prhs[6]);
    M = mxCalloc(num_materials, sizeof(Material));
    getMaterials(prhs[6], prhs[7], prhs[8], M);

    // simulation parameters
    maxScatters = (int)mxGetScalar(prhs[9]);
    n_rays = (int)mxGetScalar(prhs[10]);
    source_model = (int)mxGetScalar(prhs[11]);
    source_parameters = mxGetDoubles(prhs[12]);

    /* Set up the GSL random number generator */
    gsl_rng * my_rng = setupGSL();

    /* Put the sample and pinhole plate surface into structs */
    sample = set_up_surface(V, N, F, C, M, num_materials, ntriag_sample, sample_index);
    print_surface(&sample);

    /*
     * Create the output matrices
     * They need to be created as the transpose of what we want because of the
     * difference in indexing between MATLAB and C.
     */
    plhs[0] = mxCreateNumericMatrix(1, plate.n_detect, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, plate.n_detect*maxScatters, mxINT32_CLASS, mxREAL);

    /* Pointers to the output matrices so we may change them*/
    cntr_detected = (int32_t*)mxGetData(plhs[0]);
    numScattersRay = (int32_t*)mxGetData(plhs[2]);

    /**************************************************************************/

    /* Main implementation of the ray tracing */

    /* Loop through all the rays, tracing each one */
    // for (int i = 0; i < n_rays; i++) {
    //     Ray3D the_ray;
    //     int detected;

    //     the_ray = create_ray_source(source_parameters[0], &source_parameters[1],
    //         source_parameters[3], source_parameters[4], source_model, my_rng,
    //         source_parameters[5]);

    //     detected = trace_ray_simpleMulti(&the_ray, &killed, cntr_detected,
    //         maxScatters, Sample, Plate, the_sphere, my_rng, &detector);

    //     /*
    //      * Add the number of scattering events the ray has undergon to the
    //      * histogram. But only if it is detected.
    //      */
    //     if (detected) {
    //         int ind;
    //         ind = (detector - 1)*maxScatters + (the_ray.nScatters - 1);
    //         numScattersRay[ind]++;
    //     }
    // }

    /**************************************************************************/

    /* Output number of rays went into the detector */
    plhs[1] = mxCreateDoubleScalar(killed);

    /* Free space */
    gsl_rng_free(my_rng);
    mxFree(C);
    mxFree(M);

    return;
}


