/*
 * Copyright (c) 2020, Sam Lambrick, Dan Serment.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 */
#include "extract_inputs.h"

/*
 * Take the elements from a MATLAB cell array of strings
 * and put them in a C matrix of characters.
 * Returns how many strings were extracted.
 *
 * NB: this function will not allocate space for the character matrix.
 * You should allocate space in the same function that frees it, i.e. the caller.
 */
int get_string_cell_arr(const mxArray * cell_array, char ** strings) {
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
AnalytSphere get_sphere(const mxArray * theSphere, int index) {
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
    Material mat;
    char *names[] = {"sphere"};
    get_materials(sph_material, names, &mat);

    AnalytSphere sph = set_up_sphere(make_sphere, sphere_c, sphere_r, mat, index);
    return sph;
}


/* Extract array of material structs. Return how many were extracted. */
int get_materials(const mxArray * mat, char ** names, Material * target) {
    mxArray * field;
    if(!mxIsStruct(mat))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:material",
                          "Must be struct array");
    int n_materials = mxGetN(mat);

    for(int imat = 0; imat < n_materials; imat++) {
        field = mxGetField(mat, imat, "params");
        if(!mxIsDouble(field))
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:material",
                            "params must be doubles array.");
        int num_params = mxGetN(field);
        double * params = mxGetDoubles(field);

        field = mxGetField(mat, imat, "function");
        if(!mxIsChar(field))
            mexErrMsgIdAndTxt("MyToolbox:tracingMex:material",
                            "function must be char array.");
        char * function = mxArrayToString(field);

        target[imat] = set_up_material(names[imat], function, params, num_params);
    }
    return n_materials;
}


/**
 * Package the names, function names and parameters in the input into
 * an array of Material pointers. Return the number of materials found.
 */
int get_materials_array(const mxArray * names, const mxArray * functions,
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

    get_string_cell_arr(names, names_extr);
    get_string_cell_arr(functions, functions_extr);

    for(int idx = 0; idx < num; idx++) {
        mxArray * params_cell = mxGetCell(params, idx);
        int num_params = mxGetN(params_cell);

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

/* Extract source properties from a MATLAB array and write them to the given pointers */
void get_source(const mxArray * source, double * pinhole_r, double * pinhole_c,
                double * theta_max, double * init_angle, double * sigma) {

    const int n_params = 7;
    if(!mxIsDouble(source))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:source",
                            "Source parameters must be double array.");
    if(mxGetN(source) != n_params)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:source",
                          "Source must have exactly %d parameters.", n_params);

    double * source_parameters = mxGetDoubles(source);

    *pinhole_r = source_parameters[0];      // pinhole radius
    pinhole_c[0] = source_parameters[1];
    pinhole_c[1] = source_parameters[2];
    pinhole_c[2] = source_parameters[3];    // pinhole centre
    *theta_max = source_parameters[4];      // source max angle
    *init_angle = source_parameters[5];     // source initial angle
    *sigma = source_parameters[6];          // source std dev
}


/*
 * Extract the plate properties from MATLAB struct of length 5
 * and write to NBackWall struct.
 *
 * INPUTS:
 * - plate_opts = mxArray containing options
 * - plate_index = the surface index of the plate
 */
NBackWall get_plate(const mxArray * plate_opts, int plate_index) {
    NBackWall plate;
    mxArray * field;

    if(!mxIsStruct(plate_opts))
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate_opts",
                          "Must be struct array");
    if(mxGetN(plate_opts) != 1 || mxGetM(plate_opts) != 1)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate_opts",
                          "Must be one struct");
        
    // Get whether to represent the pinhole plate as scattering
    field = mxGetField(plate_opts, 0, "plate_represent");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:plate_opts",
                          "Representation of pinhole plate must be scalar");
    plate.plate_represent = (int)mxGetScalar(field);
    
    // Get the number of detectors
    field = mxGetField(plate_opts, 0, "n_detectors");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:plate_opts",
                          "Number of detectors must be scalar");
    plate.n_detect = (int)mxGetScalar(field);
    
    // Get the radius of the pinhole plate
    field = mxGetField(plate_opts, 0, "circle_plate_r");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:plate_opts",
                          "Radius of pinhole plate must be scalar");
    plate.circle_plate_r = mxGetScalar(field);
    
    // Get the axes of the apertures
    field = mxGetField(plate_opts, 0, "aperture_axes");
    if (!mxIsDouble(field) || mxGetN(field) != 2*plate.n_detect)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate_opts",
                          "Aperture axes must be array of 2xnumber of detectors doubles.");
    plate.aperture_axes = mxGetDoubles(field);
    
    // Get the centres of the apertures
    field = mxGetField(plate_opts, 0, "aperture_c");
    if (!mxIsDouble(field) || mxGetN(field) != 2*plate.n_detect)
        mexErrMsgIdAndTxt("MyToolbox:tracingMex:plate_opts",
                          "Aperture centres must be array of 2xnumber of detectors doubles.");
    plate.aperture_c = mxGetDoubles(field);

    // Get the number of detectors
    field = mxGetField(plate_opts, 0, "n_detectors");
    if (!mxIsScalar(field))
        mexErrMsgIdAndTxt("MuToolbox:tracingMex:plate_opts",
                          "Number of detectors must be scalar");

    plate.surf_index = plate_index;

    return plate;
}
